#include "Simulation.h"

#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscOperators.h"
#include "PetscWrappers/PetscKSP.h"
#include "PetscWrappers/PetscIS.h"

void Simulation::solveTISE()
{   
    // If we are not running the TISE, exit
    if (!getTISEStatus())
    {
        return;
    }

    // Start measuring at the very beginning
    auto total_start = MPI_Wtime();

    // Create HDF5 File to store output
    std::string filePath = getTISEOutput() + std::string("/tise.h5");
    PetscHDF5 viewer{PETSC_COMM_WORLD, filePath, FILE_MODE_WRITE};
    
    // Create Kinetic Matrix
    auto K = Matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::K,K,false);

    auto Invr2 = Matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr2,Invr2,false);
    
    auto Pot = Matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(getCtx().atom.getType(),Pot,false);

    auto S = Matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::S,S,false);

    K.AXPY(1.0, Pot, SAME_NONZERO_PATTERN);

    EPSSolver epssolver{PETSC_COMM_WORLD,getCtx().tise.getMaxIter(),getCtx().tise.getTol()};
    epssolver.setOperators(K,S);

    auto eigenvector = Vector{};
    S.setupVector(eigenvector);
    
    for (int l = 0; (l < getCtx().angular.getLmax()) && (l < getCtx().tise.getNmax()); ++l)
    {
        if (l > 0)
        {
            K.AXPY(l, Invr2,SAME_NONZERO_PATTERN);
        }
        
        PetscPrintf(PETSC_COMM_WORLD,"Solving for l = %d \n",l); 

        int reqPairs = getCtx().tise.getNmax() - l;
        epssolver.setDimensions(reqPairs);
        epssolver.solve();

        for (int convPair = 0; convPair < epssolver.getNconv(); ++convPair)
        {
            auto eigenvalue = epssolver.getEigenvalue(convPair);
            epssolver.getEigenvector(convPair,eigenvector);

            normalize(eigenvector,S);

            if (PetscRealPart(eigenvalue) > 0)
            {
                continue;
            }

            std::string eigenvectorName = std::string("psi_") + std::to_string(convPair + l + 1) + std::string("_") + std::to_string(l);
            std::string eigenvalueName = std::string("E_") + std::to_string(convPair + l + 1) + std::string("_") + std::to_string(l);

            viewer.saveValue(m_eigenvalueGroup, eigenvalueName, eigenvalue);
            viewer.saveVector(m_eigenvectorGroup,eigenvectorName, eigenvector);

            auto normVal = norm(eigenvector,S);
            PetscPrintf(PETSC_COMM_WORLD,"Eigenvector %d -> Norm(%.4f , %.4f) -> Eigenvalue(%.4f , %.4f)  \n",convPair+1,normVal.real(),normVal.imag(),eigenvalue.real(),eigenvalue.imag()); 
        }
    }
    auto total_end = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD,"Solved eigenvalue problem in:  %f seconds \n\n", total_end -  total_start);
}

void Simulation::populateAngularMatrix(AngularMatrixType Type, Matrix& matrix)
{ 
    for (int blockRow = 0; blockRow < getCtx().angular.getNblocks(); ++blockRow)
    {
        int l{};
        int m{};
        std::tie(l,m) = getCtx().angular.getBlockMap().at(blockRow);

        for (int blockCol = 0; blockCol < getCtx().angular.getNblocks(); ++blockCol)
        {
            int lprime{};
            int mprime{};
            std::tie(lprime,mprime) = getCtx().angular.getBlockMap().at(blockCol);

            switch(Type)
            {
                case AngularMatrixType::Z_INT_1:
                    if ((l == lprime + 1) && (m == mprime))
                    {
                        matrix.setValue(blockRow,blockCol,-PETSC_i * AngularElements::g(l,m));
                    }
                    else if ((l == lprime - 1)&&(m == mprime))
                    {   
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::f(l,m));
                    }
                break;
                case AngularMatrixType::Z_INT_2:
                    if ((l == lprime + 1) && (m == mprime))
                    {   
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::g(l,m) * (-l));
                    }
                    else if ((l == lprime - 1)&&(m == mprime))
                    {   
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::f(l,m) * (l+1));
                    }
                break;
                case AngularMatrixType::XY_INT_1:
                    if ((l == lprime + 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::a(l,m)/2);
                    }
                    else if ((l == lprime - 1)&&(m == mprime+1))
                    {   
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::b(l,m)/2);
                    }
                break;
                case AngularMatrixType::XY_INT_2:
                    if ((l == lprime + 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::c(l,m)/2);
                    }
                    else if ((l == lprime - 1)&&(m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::d(l,m)/2);
                    }
                break;
                case AngularMatrixType::XY_INT_3:
                    if ((l == lprime + 1) && (m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::atilde(l,m)/2);
                    }
                    else if ((l == lprime - 1)&&(m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::btilde(l,m)/2);
                    }
                break;
                case AngularMatrixType::XY_INT_4:
                    if ((l == lprime + 1) && (m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::ctilde(l,m)/2);
                    }
                    else if ((l == lprime - 1)&&(m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::dtilde(l,m)/2);
                    }
                break;
                case AngularMatrixType::X_HHG:
                    if ((l == lprime + 1) && (m == mprime - 1))
                    {   
                        matrix.setValue(blockRow,blockCol,AngularElements::charlie(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime - 1))
                    {   
                        matrix.setValue(blockRow,blockCol,AngularElements::delta(l,m));
                    }
                    else if ((l == lprime + 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol,-AngularElements::alpha(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol,-AngularElements::beta(l,m));
                    }
                break;
                case AngularMatrixType::Y_HHG:
                    if ((l == lprime + 1) && (m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::charlie(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::delta(l,m));
                    }
                    else if ((l == lprime + 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::alpha(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::beta(l,m));
                    }
                break;
                case AngularMatrixType::Z_HHG:
                    if ((l == lprime + 1) && (m == mprime))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::echo(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::foxtrot(l,m));
                    }
                break;
            }
        }
    matrix.assemble();
    }
}

void Simulation::populateRadialMatrix(RadialMatrixType Type,Matrix& matrix,bool use_ecs) 
{   
    MatrixIntegrand integrand{};

    switch(Type)
    {
        case RadialMatrixType::S:
            {
                integrand = &RadialElements::overlapIntegrand;
            }
        break;
        case RadialMatrixType::Invr2:
            {
                integrand = &RadialElements::invr2Integrand;
            }
        break;
        case RadialMatrixType::Invr:
            {
                integrand = &RadialElements::invrIntegrand;
            }
        break;
        case RadialMatrixType::Der:
            {
                integrand = &RadialElements::derIntegrand;
            }
        break;
        case RadialMatrixType::K:
            {
                integrand = &RadialElements::kineticIntegrand;
            }
        break;
        case RadialMatrixType::H:
            {
                integrand = &RadialElements::HIntegrand;
            }
        break;
    }

    if (!integrand)
    {
        return;
    }


    for (int i = matrix.getStart(); i < matrix.getEnd(); i++) 
    {
        int col_start = std::max(0, i - getCtx().basis.getOrder() + 1);
        int col_end = std::min(getCtx().basis.getNbasis(), i + getCtx().basis.getOrder());

        for (int j = col_start; j < col_end; j++) 
        {
            std::complex<double> result = getCtx().basis.integrateMatrixElement(i, j,integrand,use_ecs);
            matrix.setValue(i,j,result);
        }
    }
    matrix.assemble();
}


Matrix Simulation::kroneckerProduct(const Matrix& A, const Matrix& B,PetscInt nnz_A, PetscInt nnz_B) 
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



Vector Simulation::loadInitialState()
{   
    // Create a zero initialized vector to hold the initialte state
    auto initialState = Vector{PETSC_COMM_WORLD,PETSC_DETERMINE,getCtx().basis.getNbasis() * getCtx().angular.getNblocks()};
    initialState.setConstant(0.0);

    // Define group name and vector name to read vector from hdf5
    std::string eigenvectorGroup = "eigenvectors";
    std::string vectorName = std::string("psi_") + std::to_string(getCtx().tdse.getInitialN()) + std::string("_")  + std::to_string(getCtx().tdse.getInitialL());

    // Create hdf5 viewer for reading from file
    // Read in the vector 
    std::string filePath = getTISEOutput() + std::string("/tise.h5");
    PetscHDF5 viewer{PETSC_COMM_SELF, filePath, FILE_MODE_READ};
    auto tiseOutput = viewer.loadVector(eigenvectorGroup ,vectorName,getCtx().basis.getNbasis());


    // Extract tiseOutput to easily readable array
    const PetscScalar* tiseOutputArray;
    VecGetArrayRead(tiseOutput.get(), &tiseOutputArray);

    int blockIdx = getCtx().angular.getLMMap().at(std::make_pair(getCtx().tdse.getInitialL(),getCtx().tdse.getInitialM()));
    for (int localIdx = 0; localIdx < getCtx().basis.getNbasis(); ++localIdx)
    {
        int globalIdx = blockIdx * getCtx().basis.getNbasis() + localIdx;

        if (globalIdx >= initialState.getStart() && globalIdx < initialState.getEnd())
        {   
            initialState.setValue(globalIdx, tiseOutputArray[localIdx]);
        }
    }
    initialState.assemble();

    VecRestoreArrayRead(tiseOutput.get(), &tiseOutputArray);

    return initialState;
}


std::pair<Matrix,Matrix> Simulation::constructAtomicInteraction()
{
    Matrix I{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(), 1};
    for (int blockIdx = 0; blockIdx < getCtx().angular.getNblocks(); ++blockIdx)
    {
        I.setValue(blockIdx, blockIdx, 1.0);
    }
    I.assemble();

 
    auto totalLeft = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::S,totalLeft,true);
   
    auto totalRight = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::S,totalRight,true);

    auto K = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::K,K,true);
    
    auto Pot = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(getCtx().atom.getType(),Pot,true);

    totalLeft.AXPY(PETSC_i * getCtx().laser.getTimeSpacing() / 2.0, K, SAME_NONZERO_PATTERN);
    totalLeft.AXPY(PETSC_i * getCtx().laser.getTimeSpacing() / 2.0, Pot, SAME_NONZERO_PATTERN);

    totalRight.AXPY(-PETSC_i * getCtx().laser.getTimeSpacing() / 2.0, K, SAME_NONZERO_PATTERN);
    totalRight.AXPY(-PETSC_i * getCtx().laser.getTimeSpacing() / 2.0, Pot, SAME_NONZERO_PATTERN);


    auto firstTermLeft = kroneckerProduct(I,totalLeft,1.0,2*getCtx().basis.getDegree() + 1);
    auto firstTermRight = kroneckerProduct(I,totalRight,1.0,2*getCtx().basis.getDegree() + 1);

    Matrix modifiedI{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(), 1};
    for (int blockIdx = 0; blockIdx < getCtx().angular.getNblocks(); ++blockIdx)
    {   
        int l{};
        int m{};
        std::tie(l,m) = getCtx().angular.getBlockMap().at(blockIdx);
        modifiedI.setValue(blockIdx, blockIdx, l*(l+1) / 2.0);
    }
    modifiedI.assemble();

    auto centrifugalLeft = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr2,centrifugalLeft,true);
    
    auto centrifugalRight = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr2,centrifugalRight,true);

    centrifugalLeft *= PETSC_i * getCtx().laser.getTimeSpacing() / 2.0;
    centrifugalRight *= -PETSC_i * getCtx().laser.getTimeSpacing() / 2.0;

    auto secondTermLeft = kroneckerProduct(modifiedI, centrifugalLeft, 1.0,2*getCtx().basis.getDegree() + 1);
    auto secondTermRight = kroneckerProduct(modifiedI, centrifugalRight, 1.0,2*getCtx().basis.getDegree() + 1);

    firstTermLeft.AXPY(1.0, secondTermLeft, SAME_NONZERO_PATTERN);
    firstTermRight.AXPY(1.0, secondTermRight, SAME_NONZERO_PATTERN);

    return std::make_pair(firstTermLeft,firstTermRight);
}


Matrix Simulation::constructZInteraction()
{
    auto Der = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Der,Der,true);


    auto Invr = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr,Invr,true);

    Matrix Hlm_z_1{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),2};
    populateAngularMatrix(AngularMatrixType::Z_INT_1,Hlm_z_1);

    Matrix Hlm_z_2{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),2};
    populateAngularMatrix(AngularMatrixType::Z_INT_2,Hlm_z_2);
    
    auto ZInteraction = kroneckerProduct(Hlm_z_1, Der, 2, 2*getCtx().basis.getDegree() + 1);
    ZInteraction.AXPY(1.0,kroneckerProduct(Hlm_z_2, Invr, 2, 2*getCtx().basis.getDegree() + 1),SAME_NONZERO_PATTERN);

    return ZInteraction;
}

std::pair<Matrix,Matrix> Simulation::constructXYInteraction()
{
    auto Der = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Der,Der,true);

    auto Invr = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr,Invr,true);

    // Compute Hxy_1
    Matrix Hlm_xy_1{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),2};
    populateAngularMatrix(AngularMatrixType::XY_INT_1,Hlm_xy_1);

    Matrix Hlm_xy_2{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),2};
    populateAngularMatrix(AngularMatrixType::XY_INT_2,Hlm_xy_2);

    auto Hxy_1 = kroneckerProduct(Hlm_xy_1,Invr,2,2*getCtx().basis.getDegree() + 1);
    Hxy_1.AXPY(1.0,kroneckerProduct(Hlm_xy_2,Der,2, 2*getCtx().basis.getDegree() + 1),SAME_NONZERO_PATTERN);

    // Construct Hxy_2
    Matrix Hlm_xy_3{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),2};
    populateAngularMatrix(AngularMatrixType::XY_INT_3,Hlm_xy_3);

    Matrix Hlm_xy_4{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),2};
    populateAngularMatrix(AngularMatrixType::XY_INT_4,Hlm_xy_4);

    auto Hxy_2 = kroneckerProduct(Hlm_xy_3,Invr,2,2*getCtx().basis.getDegree() + 1);
    Hxy_2.AXPY(1.0,kroneckerProduct(Hlm_xy_4,Der,2,2*getCtx().basis.getDegree()+1),SAME_NONZERO_PATTERN);

    return std::make_pair(Hxy_1,Hxy_2);
}

Matrix Simulation::constructAtomicS()
{
    Matrix I{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(), 1};
    for (int blockIdx = 0; blockIdx < getCtx().angular.getNblocks(); ++blockIdx)
    {
        I.setValue(blockIdx, blockIdx, 1.0);
    }
    I.assemble();

    auto S = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::S,S,true);
    
    return kroneckerProduct(I,S,1,2*getCtx().basis.getDegree() + 1);
}

Matrix Simulation::constructXHHG()
{
    auto Invr2 = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr2,Invr2,false);

    Matrix Hlm_hhg_x{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),4};
    populateAngularMatrix(AngularMatrixType::X_HHG,Hlm_hhg_x);

    return kroneckerProduct(Hlm_hhg_x,Invr2,4,2*getCtx().basis.getDegree() + 1);
}

Matrix Simulation::constructYHHG()
{
    auto Invr2 = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr2,Invr2,false);

    Matrix Hlm_hhg_y{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),4};
    populateAngularMatrix(AngularMatrixType::Y_HHG,Hlm_hhg_y);

    return kroneckerProduct(Hlm_hhg_y,Invr2,4,2*getCtx().basis.getDegree() + 1);
}

Matrix Simulation::constructZHHG()
{
    auto Invr2 = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr2,Invr2,false);

    Matrix Hlm_hhg_z{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,getCtx().angular.getNblocks(),getCtx().angular.getNblocks(),4};
    populateAngularMatrix(AngularMatrixType::Z_HHG,Hlm_hhg_z);

    return kroneckerProduct(Hlm_hhg_z,Invr2,2,2*getCtx().basis.getDegree() + 1);
}

void Simulation::computeHHG(std::ofstream& hhgFile,int timeIdx,const Vector& state, const Matrix& XHHG, const Matrix& YHHG, const Matrix& ZHHG)
{
    if (getHHGStatus())
    {   
        std::complex<double> xVal{};
        std::complex<double> yVal{};
        std::complex<double> zVal{};

        if (getCtx().laser.getComponents()[0])
        {
            xVal = norm(state,XHHG);
            xVal = xVal * xVal;
        }
        if (getCtx().laser.getComponents()[1])
        {
            yVal = norm(state,YHHG);
            yVal = yVal  * yVal ;
        }
        if (getCtx().laser.getComponents()[2])
        {
            zVal = norm(state,ZHHG);
            zVal = zVal * zVal;
        }

        if (getRank()== 0)
        {
            hhgFile << timeIdx*getCtx().laser.getTimeSpacing()  << " " << xVal.real() << " "  << getCtx().laser.A(timeIdx*getCtx().laser.getTimeSpacing(),0) << " " << yVal.real() << " " << getCtx().laser.A(timeIdx*getCtx().laser.getTimeSpacing(),1)   << " " << zVal.real() << " " << getCtx().laser.A(timeIdx*getCtx().laser.getTimeSpacing(),2)  << "\n";
        }
    }
}




void Simulation::solveTDSE()
{
    if (!getTDSEStatus())
    {
        return;
    }

    // Prepate Matrices/Input state for solver
    auto start_setup = MPI_Wtime();

    auto initialState = loadInitialState();

    auto atomicS = constructAtomicS();

    auto normVal = norm(initialState,atomicS);
    PetscPrintf(PETSC_COMM_WORLD,"Initial Norm: (%.15f , %.15f) \n",normVal.real(),normVal.imag()); 

    Matrix interactionLeft{};
    Matrix interactionRight{};
    std::tie(interactionLeft,interactionRight) = constructAtomicInteraction();

    Matrix ZInteraction{};
    Matrix ZHHG{};
    if (getCtx().laser.getComponents()[2])
    {
        ZInteraction = constructZInteraction();
        if (getHHGStatus())
        {
            ZHHG = constructZHHG();
        }
    }

    Matrix Hxy_1{};
    Matrix Hxy_2{};
    Matrix XHHG{};
    Matrix YHHG{};
    if ((getCtx().laser.getComponents()[0]) || (getCtx().laser.getComponents()[1]))
    {
        std::tie(Hxy_1,Hxy_2) = constructXYInteraction();

        if (getHHGStatus() && getCtx().laser.getComponents()[0])
        {
            XHHG = constructXHHG();
        }
        if (getHHGStatus() && getCtx().laser.getComponents()[1])
        {
            YHHG = constructYHHG();
        }
    }

    PetscScalar alpha = PETSC_i * getCtx().laser.getTimeSpacing() / 2.0;

    auto end_setup = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD,"Time to setup TDSE: %f \n", end_setup-start_setup);

    // Setup HHG if necessary
    std::ofstream hhgFile;
    if (getHHGStatus() && (getRank() == 0))
    {
        hhgFile.open(std::string("misc/") + std::string("hhg_data.txt"));
        hhgFile << std::fixed << std::setprecision(15); 
    }


    // Solve TDSE
    KSPSolver ksp(PETSC_COMM_WORLD,getCtx().tdse.getMaxIter(),getCtx().tdse.getTol(),getCtx().tdse.getRestart());
    ksp.setOperators(interactionLeft);

    auto rhs = Vector{};
    interactionRight.setupVector(rhs);

    auto start_solve = MPI_Wtime();
    for (int timeIdx = 0; timeIdx < getCtx().laser.getNt(); ++timeIdx)
    {

        computeHHG(hhgFile,timeIdx,initialState,XHHG,YHHG,ZHHG);


        double tNow = timeIdx * getCtx().laser.getTimeSpacing() + getCtx().laser.getTimeSpacing() / 2.0;
        
        if (timeIdx == 0)
        {
            if (getCtx().laser.getComponents()[2])
            {
                auto Az = getCtx().laser.A(tNow,2);

                interactionLeft.AXPY(Az * alpha, ZInteraction,DIFFERENT_NONZERO_PATTERN);
                interactionRight.AXPY(-Az * alpha, ZInteraction,DIFFERENT_NONZERO_PATTERN);
            }
            if ((getCtx().laser.getComponents()[0]) || (getCtx().laser.getComponents()[1]))
            {
                auto deltaAtilde = (getCtx().laser.A(tNow,0) + PETSC_i * getCtx().laser.A(tNow,1));
                auto deltaAtildeStar = (getCtx().laser.A(tNow,0) - PETSC_i * getCtx().laser.A(tNow,1));

                interactionLeft.AXPY(alpha * deltaAtildeStar,Hxy_1,DIFFERENT_NONZERO_PATTERN);
                interactionRight.AXPY(-alpha * deltaAtildeStar,Hxy_1,DIFFERENT_NONZERO_PATTERN);

                interactionLeft.AXPY(alpha*deltaAtilde,Hxy_2,DIFFERENT_NONZERO_PATTERN);
                interactionRight.AXPY(-alpha*deltaAtilde,Hxy_2,DIFFERENT_NONZERO_PATTERN);
            }
        }
        else
        {   
            double tPrev = (timeIdx - 1) * getCtx().laser.getTimeSpacing() + getCtx().laser.getTimeSpacing() / 2.0;

            if (getCtx().laser.getComponents()[2])
            {
                auto deltaAz = (getCtx().laser.A(tNow,2) - getCtx().laser.A(tPrev,2));

                interactionLeft.AXPY( deltaAz * alpha, ZInteraction,SUBSET_NONZERO_PATTERN);
                interactionRight.AXPY( -deltaAz * alpha, ZInteraction,SUBSET_NONZERO_PATTERN);
            }
            if ((getCtx().laser.getComponents()[0]) || (getCtx().laser.getComponents()[1]))
            {
                auto deltaAtilde = ((getCtx().laser.A(tNow,0) + PETSC_i * getCtx().laser.A(tNow,1)) - (getCtx().laser.A(tPrev,0) + PETSC_i * getCtx().laser.A(tPrev,1)));
                auto deltaAtildeStar = ((getCtx().laser.A(tNow,0) - PETSC_i * getCtx().laser.A(tNow,1)) - (getCtx().laser.A(tPrev,0) - PETSC_i * getCtx().laser.A(tPrev,1)));

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

    std::string filePath = getTDSEOutput() + std::string("/tdse.h5");
    PetscHDF5 viewer(PETSC_COMM_WORLD,filePath, FILE_MODE_WRITE);
    viewer.saveVector(m_TDSEoutputGroup,m_TDSEoutputName,initialState);

    if (getHHGStatus() && (getRank() == 0))
    {
        hhgFile.close();
    }
}


void Simulation::computeBlockDistribution()
{
    if (getRank() != 0)
    {
        return;
    }
    if (!getBLOCKStatus())
    {
        return;
    }

    auto S = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::S,S,false);

    std::string filePath = getTDSEOutput() + std::string("/tdse.h5");
    PetscHDF5 viewer{PETSC_COMM_SELF,filePath,FILE_MODE_READ};
    auto finalState = viewer.loadVector(m_TDSEoutputGroup, m_TDSEoutputName,getCtx().basis.getNbasis() * getCtx().angular.getNblocks());

    std::string filename = std::string("misc") + "/block_norms.txt";

    std::ofstream outFile(filename);
    outFile << std::fixed << std::setprecision(15);

    if (!outFile)
    {
        std::cerr << "Error opening file: " << filename << '\n';
    }

    if (getCtx().observables.getProjOut())
    {
        projectOutBoundStates(finalState,S);
    }

    for (int blockIdx = 0; blockIdx < getCtx().angular.getNblocks(); ++blockIdx)
    {
        int start = blockIdx * getCtx().basis.getNbasis();

        IndexSet is{PETSC_COMM_SELF,getCtx().basis.getNbasis(), start, 1};

        Vector blockVector{};
        VecGetSubVector(finalState.get(), is.get(), &blockVector.get());

        auto normVal = innerProduct(blockVector,S,blockVector);

        VecRestoreSubVector(finalState.get(), is.get(), &blockVector.get());

        outFile << blockIdx << " " << normVal.real() << " " << normVal.imag() << '\n';
    }

    outFile.close();
}

void Simulation::projectOutBoundStates(Vector& finalState,const Matrix& S)
{   
    std::string filePath = getTISEOutput() + std::string("/tise.h5");
    PetscHDF5 viewer(PETSC_COMM_SELF,filePath,FILE_MODE_READ);

    for (int blockIdx = 0; blockIdx < getCtx().angular.getNblocks(); ++blockIdx)
    {
        int l{};
        int m{};
        std::tie(l,m) = getCtx().angular.getBlockMap().at(blockIdx);


        int start = blockIdx * getCtx().basis.getNbasis();
        
        IndexSet is{PETSC_COMM_SELF,getCtx().basis.getNbasis(),start,1};

        auto blockVector = Vector{};

        VecGetSubVector(finalState.get(),is.get(),&blockVector.get());


        for (int n = 0; n <= getCtx().tise.getNmax(); ++n)
        {
            std::string groupName = "eigenvectors";
            std::string vectorName = std::string("psi_") + std::to_string(n) + std::string("_") + std::to_string(l);

            std::string datasetName = groupName + "/" + vectorName;

            PetscBool hasDataset{};
            PetscViewerHDF5HasDataset(viewer.get(),datasetName.c_str(),&hasDataset);

            if (hasDataset)
            {
                Vector tiseState = viewer.loadVector(groupName,vectorName,getCtx().basis.getNbasis());

                PetscScalar prodVal = innerProduct(tiseState,S,blockVector);

                blockVector.AXPY(-prodVal, tiseState);
            }


        }

        VecRestoreSubVector(finalState.get(),is.get(),&blockVector.get());
    }
}

CoulombWave Simulation::computeCoulombWave(double E, int l)
{
    // Unpack necessary variables
    int Nr = getCtx().box.getNr();
    double dr = getCtx().box.getGridSpacing();

    // Initialize empty vector to store wave
    std::vector<double> wave(Nr, 0.0);

    // Compute some relevant values
    double dr2 = dr * dr;
    double k = std::sqrt(2.0 * E);
    int lterm = l * (l + 1);

    // Bootstrap the Numerov method
    wave[0] = 0.0;
    wave[1] = 1.0;

    // Numerov method to populate wave
    for (int rIdx = 2; rIdx < Nr; ++rIdx) 
    {   
        // Compute next step
        double rVal = rIdx * dr;
        double term = dr2 * (lterm/(rVal*rVal) + 2.0*getCtx().atom.potential(rVal).real() - 2.0*E);
        wave[rIdx] = wave[rIdx - 1] * (term + 2.0) - wave[rIdx - 2];

        // Check to ensure wave hasnt blown up, and normalize if it has
        if (std::abs(wave[rIdx]) > 1E10) 
        {
            double maxMag = std::abs(*std::max_element(wave.begin(), wave.end(), 
                [](double a, double b) { return std::abs(a) < std::abs(b); }));

            for (auto& value : wave)
            {
                value /= maxMag;
            }
        }
    }

    // Compute values to normalize
    double rEnd = (Nr - 2) * dr;
    double waveEnd = wave[Nr - 2];
    double dwaveEnd = (wave[Nr - 1] - wave[Nr - 3]) / (2.0 * dr);
    
    // Normalize
    double denom = k + 1.0/(k * rEnd);
    double termPsi = std::abs(waveEnd) * std::abs(waveEnd);
    double termDer = std::abs(dwaveEnd/denom) * std::abs(dwaveEnd/denom);
    double normVal = std::sqrt(termPsi + termDer);

    if (normVal > 0.0) 
    {
        for (auto& value : wave)
        {
            value /= normVal;
        }
    }

    

    // Compute values to compute phase
    std::complex<double> numerator(0.0, waveEnd);
    numerator += dwaveEnd / denom;

    const double scale = 2.0 * k * rEnd;
    const std::complex<double> denomC = std::exp(std::complex<double>(0.0, 1.0/k) * std::log(scale));
    const std::complex<double> fraction = numerator / denomC;
    const double phase = std::arg(fraction) - k * rEnd + l * M_PI/2.0;

    return CoulombWave{wave,phase};
}

std::vector<std::complex<double>> Simulation::expandState(const Vector& state)
{   
    int Nr = getCtx().box.getNr();
    double dr = getCtx().box.getGridSpacing();
    int nlm = getCtx().angular.getNblocks();
    int nbasis = getCtx().basis.getNbasis();
    int degree = getCtx().basis.getDegree();

    // Unpacked the petsc vector into a C-style array
    const std::complex<double>* stateArray;
    VecGetArrayRead(state.get(), reinterpret_cast<const PetscScalar**>(&stateArray));
    
    // Allocate vector to hold position space state
    std::vector<std::complex<double>> expanded_state(Nr * nlm);


    // Evaluate all bspline basis functions in their nonzero intervals
    for (int bsplineIdx = 0; bsplineIdx < nbasis; ++bsplineIdx)
    {   
        // Get the start and end of the interval where the bspline is nonzero
        std::complex<double> start = getCtx().basis.getKnots()[bsplineIdx];
        std::complex<double> end = getCtx().basis.getKnots()[bsplineIdx+degree+1];

        // Initialize vectors to store nonzero bspline values and the position indices where they occur
        std::vector<std::complex<double>> bsplineEval{};
        std::vector<int> bsplineEvalIndices{};

        // Loop over all grid points
        for (int rIdx = 0; rIdx < Nr; ++rIdx)
        {   
            // Compute position
            double r = rIdx*dr;

            // See if position is in nonzero interval for this spline, if so evaluate
            if (r >= start.real() && r < end.real())
            {
                std::complex<double> bsplineVal = BSplines::B(bsplineIdx,degree,r,getCtx().basis.getKnots());
                bsplineEval.push_back(bsplineVal);
                bsplineEvalIndices.push_back(rIdx);
            }
        }

        // Loop over each block 
        for (int blockIdx = 0; blockIdx < nlm; ++blockIdx)
        {   
            // Get the bspline coefficients
            std::complex<double> coeff = stateArray[blockIdx*nbasis + bsplineIdx];

            // Loop over all grid points and add contribution to the expanded state for this block
            for (size_t rSubIdx = 0; rSubIdx < bsplineEval.size(); ++rSubIdx)
            {
                expanded_state[blockIdx*Nr + bsplineEvalIndices[rSubIdx]] += coeff*bsplineEval[rSubIdx];
            }
        }
    }
    return  expanded_state;
}

std::pair<std::map<lmPair,std::vector<std::complex<double>>>,std::map<std::pair<double, int>,double>> Simulation::computePartialSpectra(const std::vector<std::complex<double>>& expanded_state)
{
    std::map<lmPair,std::vector<std::complex<double>>> partialSpectra;

    std::map<std::pair<double, int>,double> phases;

    int nlm = getCtx().angular.getNblocks();
    int Ne = getCtx().observables.getNe();
    double Emin = getCtx().observables.getEmin();
    int Nr = getCtx().box.getNr();
    double dr = getCtx().box.getGridSpacing();

    // Allocate space for each partial spectra
    for (int blockIdx = 0; blockIdx < nlm; ++blockIdx)
    {
        std::pair<int,int> lmPair = getCtx().angular.getBlockMap().at(blockIdx);
        int l = lmPair.first;
        int m = lmPair.second;
        partialSpectra[std::make_pair(l, m)].reserve(Ne); 
    }


    for (int EIdx = 1; EIdx <= Ne; ++EIdx)
    {
        for (int blockIdx = 0; blockIdx < nlm; ++blockIdx)
        {
            std::pair<int,int> lm_pair = getCtx().angular.getBlockMap().at(blockIdx);
            int l = lm_pair.first;
            int m = lm_pair.second;
            CoulombWave coulombResult = computeCoulombWave(EIdx*Emin, l);
            
            phases[std::make_pair(EIdx*Emin,l)] = coulombResult.phase;


            auto start = expanded_state.begin() + Nr*blockIdx;  
            auto end = expanded_state.begin() + Nr*(blockIdx+1);    

            // Extract subvector
            std::vector<std::complex<double>> blockVector(start, end);  

            // We want to form the position space inner product aka integrate. 
            // So first we to a pointwise mult. Since coulomb wave is real no need for complex conjugation

            for (int rIdx = 0; rIdx < Nr; ++rIdx)
            {
                blockVector[rIdx] = blockVector[rIdx] * coulombResult.wave[rIdx];
            }

            std::complex<double> I = SimpsonsMethod(blockVector,dr);  
            partialSpectra[std::make_pair(l,m)].push_back(I);
        }
    }
    return std::make_pair(partialSpectra,phases);
}

void Simulation::computeAngleIntegrated(const std::map<lmPair,std::vector<std::complex<double>>>& partialSpectra)
{   
    int Ne = getCtx().observables.getNe();
    int nlm = getCtx().angular.getNblocks();
    double Emin = getCtx().observables.getEmin();

    std::ofstream pesFiles("misc/pes.txt");
    std::vector<std::complex<double>> pes(Ne);
    
    for (int blockIdx = 0; blockIdx < nlm; ++blockIdx)
    {   
        lmPair lm_pair = getCtx().angular.getBlockMap().at(blockIdx);
        int l = lm_pair.first;
        int m = lm_pair.second;

        std::vector<std::complex<double>> partialSpectrum = partialSpectra.at(std::make_pair(l,m));

        std::vector<std::complex<double>> magSq(Ne);
        for (size_t idx = 0; idx < partialSpectrum.size(); ++idx)
        {
            magSq[idx] = partialSpectrum[idx] * std::conj(partialSpectrum[idx]);
        }

        for (size_t idx = 0; idx < magSq.size(); ++idx)
        {
            pes[idx] += magSq[idx];
        }
    }
    

    for (size_t idx = 1; idx < pes.size(); ++idx)
    {   
        std::complex<double> val = pes[idx];
        val /= ((2*M_PI)*(2*M_PI)*(2*M_PI));
        pesFiles << idx*Emin << " " << val.real() << " " << "\n";
    }

    pesFiles.close();
}

void Simulation::computeAngleResolved(const std::map<lmPair,std::vector<std::complex<double>>>& partialSpectra,std::map<std::pair<double, int>,double> phases)
{   
    int Ne = getCtx().observables.getNe();
    double Emin = getCtx().observables.getEmin();

        std::ofstream padFiles("misc/pad.txt");
        std::vector<double> theta_range;
        std::vector<double> phi_range;

        if (getCtx().observables.getSlice() == "XZ")
        {
            for (double theta = 0; theta <= M_PI; theta += 0.01) 
            {
                theta_range.push_back(theta);
            }

            phi_range  = {0.0,M_PI};
        }
        if (getCtx().observables.getSlice() == "XY")
        {
            theta_range = {M_PI/ 2.0};

            for (double phi = 0; phi < 2.0*M_PI; phi += 0.01) 
            {
                phi_range.push_back(phi);
            }

        }
        for (int EIdx = 1; EIdx <= Ne; ++EIdx)
        {   

            double E = EIdx*Emin;
            double k = std::sqrt(2.0*E);

            for (auto& theta : theta_range)
            {
                for (auto& phi : phi_range)
                {   
                    std::complex<double> pad_amplitude {};
                    for (const auto& pair: partialSpectra)
                    {
                        int l = pair.first.first;
                        int m = pair.first.second;

                        std::complex<double> sph_term = compute_Ylm(l,m,theta,phi);

                        double partial_phase = phases[std::make_pair(E,l)];
                        std::complex<double> partial_amplitude = pair.second[EIdx - 1];

                        std::complex<double> phase_factor = std::exp(std::complex<double>(0.0,3*l*M_PI/2.0 + partial_phase));

                        pad_amplitude += sph_term*phase_factor*partial_amplitude;

                    }

                    double pad_prob = std::norm(pad_amplitude);
                    pad_prob/=((2*M_PI)*(2*M_PI)*(2*M_PI));
                    pad_prob/=k;

                    padFiles << E << " " << theta << " " << phi << " " << pad_prob << "\n";
                }
            }
        }
}

void Simulation::computePhotoelectronSpectrum()
{
    if (getRank() != 0)
    {
        return;
    }

    if (!getPESStatus())
    {
        return;
    }

    auto start_pes = MPI_Wtime();

    auto S = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::S,S,false);

    std::string filePath = getTDSEOutput() + std::string("/tdse.h5");
    PetscHDF5 viewer{PETSC_COMM_SELF,filePath,FILE_MODE_READ};
    std::string groupName = "";
    std::string vectorName = "psiFinal";
    auto finalState = viewer.loadVector(groupName, vectorName,getCtx().basis.getNbasis() * getCtx().angular.getNblocks());

    projectOutBoundStates(finalState,S);

    std::vector<std::complex<double>> expandedState = expandState(finalState);


    std::map<lmPair,std::vector<std::complex<double>>> partialSpectra;
    std::map<std::pair<double, int>,double> phases;

    std::tie(partialSpectra,phases) = computePartialSpectra(expandedState);

    computeAngleIntegrated(partialSpectra);
    computeAngleResolved(partialSpectra,phases);

    auto end_pes = MPI_Wtime();
    PetscPrintf(PETSC_COMM_SELF,"Computed Photoelectron Spectrum:  %f seconds \n\n", end_pes -  start_pes);
}

double Simulation::computeBoundPopulation(int n_bound, int l_bound, const Vector& state)
{

    if ((l_bound >= n_bound) || (n_bound > getCtx().tise.getNmax()))
    {
        return 0.0;
    }
    

    double probability = 0;

    std::string filePath = getTISEOutput() + std::string("/tise.h5;");
    PetscHDF5 viewer{PETSC_COMM_SELF, filePath, FILE_MODE_READ};

    std::string eigenvectorGroup = "eigenvectors";
    std::string vectorName = std::string("psi_") + std::to_string(n_bound) + std::string("_")  + std::to_string(l_bound);
    auto boundState = viewer.loadVector(eigenvectorGroup, vectorName, getCtx().basis.getNbasis());

  


   

    auto S = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::S,S,false);

    for (int blockIdx = 0; blockIdx < getCtx().angular.getNblocks(); ++blockIdx)
    {


        std::pair<int, int> lmPair = getCtx().angular.getBlockMap().at(blockIdx);
        int l = lmPair.first;
        

        if (!(l == l_bound))
        {   
            continue;
        }

       

        int start = blockIdx * getCtx().basis.getNbasis();
        auto is = IndexSet{PETSC_COMM_SELF, getCtx().basis.getNbasis(), start, 1};

        Vector stateBlock{};
        VecGetSubVector(state.get(), is.get(), &stateBlock.get());


        std::complex<double> projection = innerProduct(boundState,S,stateBlock);
        probability += std::norm(projection);

        VecRestoreSubVector(state.get(), is.get(), &stateBlock.get());
    }
    return probability;
}

void Simulation::computeBoundDistribution()
{
    if ((getRank() != 0) || (!getBOUNDStatus()))
    {
        return;
    }

    std::ofstream file(std::string("misc/") + std::string("bound_pops.txt"));
    file << std::fixed << std::setprecision(15);

    std::string filePath = getTDSEOutput() + std::string("/tdse.h5");
    PetscHDF5 viewer{PETSC_COMM_SELF,filePath,FILE_MODE_READ};
    std::string groupName = "";
    std::string vectorName = "psiFinal";
    auto finalState = viewer.loadVector(groupName, vectorName,getCtx().basis.getNbasis() * getCtx().angular.getNblocks());

    for (int n = 1; n <= getCtx().tise.getNmax(); ++n)
    {
        for (int l = 0; l < n; ++l)
        {   
            double pop = computeBoundPopulation(n,l,finalState);
            file << n << " " << l << " " << pop << '\n';
        }
    }
    file.close();
}



// void Observables::printConfiguration(int rank)
// {
//     if (rank == 0)
//     {
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
//         std::cout << "Observables Configuration: " << "\n\n";
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
//         std::cout << "Block :  projOutBound: " << getProjOut() << "\n\n";
//     }
// }
