#include "TDSE.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscMat.h"
#include "PetscWrappers/PetscOperators.h"
#include "PetscWrappers/PetscKSP.h"

Matrix TDSE::kroneckerProduct(const Matrix& A, const Matrix& B,PetscInt nnz_A, PetscInt nnz_B) 
{
    

    // Get matrix dimensions
    PetscInt am, an, bm, bn;
    MatGetSize(A.get(), &am, &an); 
    MatGetSize(B.get(), &bm, &bn); 

    // Compute dimensions of Kronecker product matrix
    PetscInt cm = am * bm;
    PetscInt cn = an * bn;

    // Access internal csr data of A and B
    const PetscInt *ai, *aj, *bi, *bj;
    const PetscScalar *aa, *ba;
    MatGetRowIJ(A.get(), 0, PETSC_FALSE, PETSC_FALSE, &am, &ai, &aj, NULL);
    MatGetRowIJ(B.get(), 0, PETSC_FALSE, PETSC_FALSE, &bm, &bi, &bj, NULL);
    MatSeqAIJGetArrayRead(A.get(), &aa); 
    MatSeqAIJGetArrayRead(B.get(), &ba); 

    // Estimate upper bound for nnz_per row
    PetscInt nnz_C = nnz_A * nnz_B;


    // Create parallel matrix C
    auto C = Matrix{PETSC_COMM_WORLD,PETSC_DETERMINE,PETSC_DETERMINE,cm,cn,nnz_C};

    // Preallocate matrix C for local rows
    PetscInt local_nnz = 0;
    PetscInt *ci = new PetscInt[C.getEnd() - C.getStart() + 1];
    PetscInt *cj = nullptr;
    PetscScalar *cv = nullptr;


    ci[0] = 0;
    for (PetscInt iC = C.getStart(); iC < C.getEnd(); ++iC) {
        PetscInt iA = iC / bm; // Row index in A
        PetscInt iB = iC % bm; // Row index in B


        for (PetscInt n = ai[iA]; n < ai[iA + 1]; ++n) {
            for (PetscInt q = bi[iB]; q < bi[iB + 1]; ++q) {
                local_nnz++;
            }
        }


        ci[iC - C.getStart() + 1] = local_nnz;
    }


    cj = new PetscInt[local_nnz];
    cv = new PetscScalar[local_nnz];


    local_nnz = 0;
    for (PetscInt iC = C.getStart(); iC < C.getEnd(); ++iC) {
        PetscInt iA = iC / bm;
        PetscInt iB = iC % bm;


        for (PetscInt n = ai[iA]; n < ai[iA + 1]; ++n) {
            PetscInt colA = aj[n];
            PetscScalar valA = aa[n];


            for (PetscInt q = bi[iB]; q < bi[iB + 1]; ++q) {
                cj[local_nnz] = colA * bn + bj[q];
                cv[local_nnz] = valA * ba[q];
                local_nnz++;
            }
        }
    }

    

    MatMPIAIJSetPreallocationCSR(C.get(), ci, cj, cv); 


    // Clean up
    delete[] ci;
    delete[] cj;
    delete[] cv;


    MatRestoreRowIJ(A.get(), 0, PETSC_FALSE, PETSC_FALSE, &am, &ai, &aj,nullptr); 
    MatRestoreRowIJ(B.get(), 0, PETSC_FALSE, PETSC_FALSE, &bm, &bi, &bj,nullptr); 
    MatSeqAIJRestoreArrayRead(A.get(), &aa); 
    MatSeqAIJRestoreArrayRead(B.get(), &ba); 


    // Assemble the matrix
    C.assemble();

    return C;
}

void TDSE::printConfiguration(int rank)
{
    if (rank == 0)
    {
        std::cout << "TDSE Configuration: " << "\n\n";
        std::cout << "status: " << getStatus() << "\n\n";
        std::cout << "maxIter: " << getMaxIter() << "\n\n";
        std::cout << "outputPath: " << getOutputPath() << "\n\n";
        std::cout << "Restart: " << getRestart() << "\n\n";
        std::cout << "Initial nlm: " << "n = " << getInitialNLM()[0]  << " l = " << getInitialNLM()[1]  << " m = " << getInitialNLM()[2] << "\n\n";

    }
}


Vector TDSE::loadInitialState(const TISE& tise,const BSpline& bspline, const Angular& angular)
{   
    // Create a zero initialized vector to hold the initialte state
    auto initialState = Vector{PETSC_COMM_WORLD,PETSC_DETERMINE,bspline.getNbasis() * angular.getNlm()};
    initialState.setConstant(0.0);

    // Define group name and vector name to read vector from hdf5
    std::string eigenvectorGroup = "eigenvectors";
    std::string vectorName = std::string("psi_") + std::to_string(getInitialNLM()[0]) + std::string("_")  + std::to_string(getInitialNLM()[1]);

    // Create hdf5 viewer for reading from file
    // Read in the vector 
    PetscHDF5 viewer{PETSC_COMM_SELF, tise.getOutputPath(), FILE_MODE_READ};
    auto tiseOutput = viewer.loadVector(eigenvectorGroup ,vectorName,bspline.getNbasis());


    // Extract tiseOutput to easily readable array
    const PetscScalar* tiseOutputArray;
    VecGetArrayRead(tiseOutput.get(), &tiseOutputArray);

    int blockIdx = angular.getLMMap().at(std::make_pair(getInitialNLM()[1],getInitialNLM()[2]));
    for (int localIdx = 0; localIdx < bspline.getNbasis(); ++localIdx)
    {
        int globalIdx = blockIdx * bspline.getNbasis() + localIdx;

        if (globalIdx >= initialState.getStart() && globalIdx < initialState.getEnd())
        {   
            initialState.setValue(globalIdx, tiseOutputArray[localIdx]);
        }
    }
    initialState.assemble();

    VecRestoreArrayRead(tiseOutput.get(), &tiseOutputArray);

    return initialState;
}

std::pair<Matrix,Matrix> TDSE::constructAtomicInteraction(const BSpline& bspline, const Angular& angular, const Atom& atom, const Laser& laser)
{
    Matrix I{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(), 1};
    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {
        I.setValue(blockIdx, blockIdx, 1.0);
    }
    I.assemble();

    auto totalLeft = bspline.PopulateMatrix(PETSC_COMM_SELF,&BSpline::overlapIntegrand, true);
    auto totalRight = bspline.PopulateMatrix(PETSC_COMM_SELF,&BSpline::overlapIntegrand, true);

    auto K = bspline.PopulateMatrix(PETSC_COMM_SELF,&BSpline::kineticIntegrand, true);
    Matrix Pot{};
    if (atom.getSpecies() == "H")
    {
        Pot = bspline.PopulateMatrix(PETSC_COMM_SELF,&BSpline::HIntegrand, true);
    }

    totalLeft.AXPY(PETSC_i * laser.getTimeSpacing() / 2.0, K, SAME_NONZERO_PATTERN);
    totalLeft.AXPY(PETSC_i * laser.getTimeSpacing() / 2.0, Pot, SAME_NONZERO_PATTERN);

    totalRight.AXPY(-PETSC_i * laser.getTimeSpacing() / 2.0, K, SAME_NONZERO_PATTERN);
    totalRight.AXPY(-PETSC_i * laser.getTimeSpacing() / 2.0, Pot, SAME_NONZERO_PATTERN);


    auto firstTermLeft = kroneckerProduct(I,totalLeft,1.0,2*bspline.getDegree() + 1);
    auto firstTermRight = kroneckerProduct(I,totalRight,1.0,2*bspline.getDegree() + 1);

    Matrix modifiedI{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(), 1};
    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {   
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockIdx);
        modifiedI.setValue(blockIdx, blockIdx, l*(l+1) / 2.0);
    }
    modifiedI.assemble();

    auto centrifugalLeft = bspline.PopulateMatrix(PETSC_COMM_SELF,&BSpline::invr2Integrand, true);
    auto centrifugalRight = bspline.PopulateMatrix(PETSC_COMM_SELF,&BSpline::invr2Integrand, true);

    centrifugalLeft *= PETSC_i * laser.getTimeSpacing() / 2.0;
    centrifugalRight *= -PETSC_i * laser.getTimeSpacing() / 2.0;

    auto secondTermLeft = kroneckerProduct(modifiedI, centrifugalLeft, 1.0,2*bspline.getDegree() + 1);
    auto secondTermRight = kroneckerProduct(modifiedI, centrifugalRight, 1.0,2*bspline.getDegree() + 1);

    firstTermLeft.AXPY(1.0, secondTermLeft, SAME_NONZERO_PATTERN);
    firstTermRight.AXPY(1.0, secondTermRight, SAME_NONZERO_PATTERN);

    return std::make_pair(firstTermLeft,firstTermRight);
}

Matrix TDSE::constructZInteraction(const BSpline& bspline, const Angular& angular)
{
    auto Invr = bspline.PopulateMatrix(PETSC_COMM_SELF, &BSpline::invrIntegrand, true);
    auto Der = bspline.PopulateMatrix(PETSC_COMM_SELF, &BSpline::derIntegrand, true);

    Matrix Hlm_z_1{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
    for (int blockRow = 0; blockRow < angular.getNlm(); ++blockRow)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockRow);

        for (int blockCol = 0; blockCol < angular.getNlm(); ++blockCol)
        {
            int lprime{};
            int mprime{};
            std::tie(lprime,mprime) = angular.getBlockMap().at(blockCol);


            if ((l == lprime + 1) && (m == mprime))
            {
                Hlm_z_1.setValue(blockRow,blockCol, -PETSC_i * Angular::Coupling::g(l,m));
            }
            if ((l == lprime - 1) && (m == mprime))
            {
                Hlm_z_1.setValue(blockRow,blockCol, -PETSC_i * Angular::Coupling::f(l,m));
            }
        }
    }
    Hlm_z_1.assemble();

    Matrix Hlm_z_2{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
    for (int blockRow = 0; blockRow < angular.getNlm(); ++blockRow)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockRow);

        for (int blockCol = 0; blockCol < angular.getNlm(); ++blockCol)
        {
            int lprime{};
            int mprime{};
            std::tie(lprime,mprime) = angular.getBlockMap().at(blockCol);


            if ((l == lprime + 1) && (m == mprime))
            {
                Hlm_z_2.setValue(blockRow,blockCol, -PETSC_i * Angular::Coupling::g(l,m) * (-l));
            }
            if ((l == lprime - 1) && (m == mprime))
            {
                Hlm_z_2.setValue(blockRow,blockCol, -PETSC_i * Angular::Coupling::f(l,m) * (l+1));
            }
        }
    }
    Hlm_z_2.assemble();

    auto ZInteraction = kroneckerProduct(Hlm_z_1, Der, 2, 2*bspline.getDegree() + 1);
    ZInteraction.AXPY(1.0,kroneckerProduct(Hlm_z_2, Invr, 2, 2*bspline.getDegree() + 1),SAME_NONZERO_PATTERN);

    return ZInteraction;
}

std::pair<Matrix,Matrix> TDSE::constructXYInteraction(const BSpline& bspline, const Angular& angular)
{
    auto Invr = bspline.PopulateMatrix(PETSC_COMM_SELF, &BSpline::invrIntegrand, true);
    auto Der = bspline.PopulateMatrix(PETSC_COMM_SELF, &BSpline::derIntegrand, true);



    // Compute Hxy_1
    Matrix Hlm_xy_1{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
    for (int blockRow = 0; blockRow < angular.getNlm(); ++blockRow)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockRow);

        for (int blockCol = 0; blockCol < angular.getNlm(); ++blockCol)
        {
            int lprime{};
            int mprime{};
            std::tie(lprime,mprime) = angular.getBlockMap().at(blockCol);


            if ((l == lprime + 1) && (m == mprime + 1))
            {
                Hlm_xy_1.setValue(blockRow,blockCol, PETSC_i * Angular::Coupling::a(l,m) / 2.0);
            }
            if ((l == lprime - 1) && (m == mprime + 1))
            {
                Hlm_xy_1.setValue(blockRow,blockCol, PETSC_i * Angular::Coupling::b(l,m) / 2.0);
            }
        }
    }
    Hlm_xy_1.assemble();

    Matrix Hlm_xy_2{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
    for (int blockRow = 0; blockRow < angular.getNlm(); ++blockRow)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockRow);

        for (int blockCol = 0; blockCol < angular.getNlm(); ++blockCol)
        {
            int lprime{};
            int mprime{};
            std::tie(lprime,mprime) = angular.getBlockMap().at(blockCol);


            if ((l == lprime + 1) && (m == mprime + 1))
            {
                Hlm_xy_2.setValue(blockRow,blockCol, PETSC_i * Angular::Coupling::c(l,m) / 2.0);
            }
            if ((l == lprime - 1) && (m == mprime + 1))
            {
                Hlm_xy_2.setValue(blockRow,blockCol, -PETSC_i * Angular::Coupling::d(l,m) / 2.0);
            }
        }
    }
    Hlm_xy_2.assemble();

    auto Hxy_1 = kroneckerProduct(Hlm_xy_1,Invr,2,2*bspline.getDegree() + 1);
    Hxy_1.AXPY(1.0,kroneckerProduct(Hlm_xy_2,Der,2, 2*bspline.getDegree() + 1),SAME_NONZERO_PATTERN);

    // Construct Hxy_2
    Matrix Hlm_xy_3{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
    for (int blockRow = 0; blockRow < angular.getNlm(); ++blockRow)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockRow);

        for (int blockCol = 0; blockCol < angular.getNlm(); ++blockCol)
        {
            int lprime{};
            int mprime{};
            std::tie(lprime,mprime) = angular.getBlockMap().at(blockCol);


            if ((l == lprime + 1) && (m == mprime - 1))
            {
                Hlm_xy_3.setValue(blockRow,blockCol, PETSC_i * Angular::Coupling::atilde(l,m) / 2.0);
            }
            if ((l == lprime - 1) && (m == mprime - 1))
            {
                Hlm_xy_3.setValue(blockRow,blockCol, PETSC_i * Angular::Coupling::btilde(l,m) / 2.0);
            }
        }
    }
    Hlm_xy_3.assemble();

    Matrix Hlm_xy_4{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
    for (int blockRow = 0; blockRow < angular.getNlm(); ++blockRow)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockRow);

        for (int blockCol = 0; blockCol < angular.getNlm(); ++blockCol)
        {
            int lprime{};
            int mprime{};
            std::tie(lprime,mprime) = angular.getBlockMap().at(blockCol);


            if ((l == lprime + 1) && (m == mprime - 1))
            {
                Hlm_xy_4.setValue(blockRow,blockCol, -PETSC_i * Angular::Coupling::ctilde(l,m) / 2.0);
            }
            if ((l == lprime - 1) && (m == mprime - 1))
            {
                Hlm_xy_4.setValue(blockRow,blockCol, PETSC_i * Angular::Coupling::dtilde(l,m) / 2.0);
            }
        }
    }
    Hlm_xy_4.assemble();

    auto Hxy_2 = kroneckerProduct(Hlm_xy_3,Invr,2,2*bspline.getDegree() + 1);
    Hxy_2.AXPY(1.0,kroneckerProduct(Hlm_xy_4,Der,2,2*bspline.getDegree()+1),SAME_NONZERO_PATTERN);

    return std::make_pair(Hxy_1,Hxy_2);
}

Matrix TDSE::constructAtomicS(const BSpline& bspline, const Angular& angular)
{
    Matrix I{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(), 1};
    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {
        I.setValue(blockIdx, blockIdx, 1.0);
    }
    I.assemble();

    auto S = bspline.PopulateMatrix(PETSC_COMM_SELF,&BSpline::overlapIntegrand, true);
    
    return kroneckerProduct(I,S,1,2*bspline.getDegree() + 1);
}


void TDSE::solve(const TISE& tise,const BSpline& bspline, const Angular& angular, const Atom& atom, const Laser& laser)
{
    if (!getStatus())
    {
        return;
    }

    // Prepate Matrices/Input state for solver
    auto start_setup = MPI_Wtime();

    auto initialState = loadInitialState(tise,bspline,angular);

    auto atomicS = constructAtomicS(bspline,angular);

    auto normVal = norm(initialState,atomicS);
    PetscPrintf(PETSC_COMM_WORLD,"Initial Norm: (%.15f , %.15f) \n",normVal.real(),normVal.imag()); 

    Matrix interactionLeft{};
    Matrix interactionRight{};
    std::tie(interactionLeft,interactionRight) = constructAtomicInteraction(bspline,angular,atom,laser);

    Matrix ZInteraction{};
    if (laser.getComponents()[2])
    {
        ZInteraction = constructZInteraction(bspline,angular);
    }

    Matrix Hxy_1{};
    Matrix Hxy_2{};
    if ((laser.getComponents()[0]) || (laser.getComponents()[1]))
    {
        std::tie(Hxy_1,Hxy_2) = constructXYInteraction(bspline,angular);
    }

    PetscScalar alpha = PETSC_i * laser.getTimeSpacing() / 2.0;

    auto end_setup = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD,"Time to setup TDSE: %f \n", end_setup-start_setup);


    // Solve TDSE
    KSPSolver ksp(PETSC_COMM_WORLD,getMaxIter(),getTol(),getRestart());
    ksp.setOperators(interactionLeft);

    auto rhs = Vector{};
    interactionRight.setupVector(rhs);

    auto start_solve = MPI_Wtime();
    for (int timeIdx = 0; timeIdx < laser.getNt(); ++timeIdx)
    {
        double tNow = timeIdx * laser.getTimeSpacing() + laser.getTimeSpacing() / 2.0;
        
        if (timeIdx == 0)
        {
            if (laser.getComponents()[2])
            {
                auto Az = laser.A(tNow,2);;

                interactionLeft.AXPY(Az * alpha, ZInteraction,DIFFERENT_NONZERO_PATTERN);
                interactionRight.AXPY(-Az * alpha, ZInteraction,DIFFERENT_NONZERO_PATTERN);
            }
            if ((laser.getComponents()[0]) || (laser.getComponents()[1]))
            {
                auto deltaAtilde = (laser.A(tNow,0) + PETSC_i * laser.A(tNow,1));
                auto deltaAtildeStar = (laser.A(tNow,0) - PETSC_i * laser.A(tNow,1));;

                interactionLeft.AXPY(alpha * deltaAtildeStar,Hxy_1,DIFFERENT_NONZERO_PATTERN);
                interactionRight.AXPY(-alpha * deltaAtildeStar,Hxy_1,DIFFERENT_NONZERO_PATTERN);

                interactionLeft.AXPY(alpha*deltaAtilde,Hxy_2,DIFFERENT_NONZERO_PATTERN);
                interactionRight.AXPY(-alpha*deltaAtilde,Hxy_2,DIFFERENT_NONZERO_PATTERN);
            }
        }
        else
        {   
            double tPrev = (timeIdx - 1) * laser.getTimeSpacing() + laser.getTimeSpacing() / 2.0;

            if (laser.getComponents()[2])
            {
                auto deltaAz = (laser.A(tNow,2) - laser.A(tPrev,2));

                interactionLeft.AXPY( deltaAz * alpha, ZInteraction,SUBSET_NONZERO_PATTERN);
                interactionRight.AXPY( -deltaAz * alpha, ZInteraction,SUBSET_NONZERO_PATTERN);
            }
            if ((laser.getComponents()[0]) || (laser.getComponents()[1]))
            {
                auto deltaAtilde = ((laser.A(tNow,0) + PETSC_i * laser.A(tNow,1)) - (laser.A(tPrev,0) + PETSC_i * laser.A(tPrev,1)));
                auto deltaAtildeStar = ((laser.A(tNow,0) - PETSC_i * laser.A(tNow,1)) - (laser.A(tPrev,0) - PETSC_i * laser.A(tPrev,1)));

                interactionLeft.AXPY(alpha * deltaAtildeStar,Hxy_1,SUBSET_NONZERO_PATTERN);
                interactionRight.AXPY(-alpha * deltaAtildeStar,Hxy_1,SUBSET_NONZERO_PATTERN);

                interactionLeft.AXPY(alpha*deltaAtilde,Hxy_2,SUBSET_NONZERO_PATTERN);
                interactionRight.AXPY(-alpha*deltaAtilde,Hxy_2,SUBSET_NONZERO_PATTERN);
            }
        }


        interactionRight.matMult(initialState,rhs);
        ksp.solve(rhs,initialState);
        
    }

    
    auto end_solve = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD,"Time to solve TDSE: %f \n", end_solve-start_solve);

    normVal = norm(initialState,atomicS);
    PetscPrintf(PETSC_COMM_WORLD,"Final Norm: (%.15f , %.15f) \n",normVal.real(),normVal.imag()); 

    std::string outputGroup = "";
    std::string outputName = "psiFinal";

    PetscHDF5 viewer(PETSC_COMM_WORLD,getOutputPath(), FILE_MODE_WRITE);
    viewer.saveVector(outputGroup,outputName,initialState);

    
}
