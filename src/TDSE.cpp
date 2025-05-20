#include "TDSE.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscMat.h"
#include "PetscWrappers/PetscOperators.h"
#include "PetscWrappers/PetscKSP.h"
#include "MatrixElements.h"

#include "MatrixElements.h"
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

// void TDSE::printConfiguration(int rank)
// {
//     if (rank == 0)
//     {
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
//         std::cout << "TDSE Configuration: " << "\n\n";
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
//         std::cout << "status: " << getStatus() << "\n\n";
//         std::cout << "maxIter: " << getMaxIter() << "\n\n";
//         std::cout << "outputPath: " << getOutputPath() << "\n\n";
//         std::cout << "Restart: " << getRestart() << "\n\n";
//         std::cout << "Initial nlm: " << "n = " << getInitialNLM()[0]  << " l = " << getInitialNLM()[1]  << " m = " << getInitialNLM()[2] << "\n\n";

//     }
//}

// Vector TDSE::loadInitialState(const TISE& tise,const Basis& basis, const Angular& angular)
// {   
//     // Create a zero initialized vector to hold the initialte state
//     auto initialState = Vector{PETSC_COMM_WORLD,PETSC_DETERMINE,basis.getNbasis() * angular.getNlm()};
//     initialState.setConstant(0.0);

//     // Define group name and vector name to read vector from hdf5
//     std::string eigenvectorGroup = "eigenvectors";
//     std::string vectorName = std::string("psi_") + std::to_string(getInitialNLM()[0]) + std::string("_")  + std::to_string(getInitialNLM()[1]);

//     // Create hdf5 viewer for reading from file
//     // Read in the vector 
//     PetscHDF5 viewer{PETSC_COMM_SELF, tise.getOutputPath(), FILE_MODE_READ};
//     auto tiseOutput = viewer.loadVector(eigenvectorGroup ,vectorName,basis.getNbasis());


//     // Extract tiseOutput to easily readable array
//     const PetscScalar* tiseOutputArray;
//     VecGetArrayRead(tiseOutput.get(), &tiseOutputArray);

//     int blockIdx = angular.getLMMap().at(std::make_pair(getInitialNLM()[1],getInitialNLM()[2]));
//     for (int localIdx = 0; localIdx < basis.getNbasis(); ++localIdx)
//     {
//         int globalIdx = blockIdx * basis.getNbasis() + localIdx;

//         if (globalIdx >= initialState.getStart() && globalIdx < initialState.getEnd())
//         {   
//             initialState.setValue(globalIdx, tiseOutputArray[localIdx]);
//         }
//     }
//     initialState.assemble();

//     VecRestoreArrayRead(tiseOutput.get(), &tiseOutputArray);

//     return initialState;
// }

// std::pair<Matrix,Matrix> TDSE::constructAtomicInteraction(const Basis& basis, const Angular& angular, const Atom& atom, const Laser& laser)
// {
//     Matrix I{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(), 1};
//     for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
//     {
//         I.setValue(blockIdx, blockIdx, 1.0);
//     }
//     I.assemble();

 
//     auto totalLeft = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::S,totalLeft,basis,true);
   
//     auto totalRight = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::S,totalRight,basis,true);

//     auto K = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::K,K,basis,true);
    
//     auto Pot = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(atom.getType(),Pot,basis,true);

//     totalLeft.AXPY(PETSC_i * laser.getTimeSpacing() / 2.0, K, SAME_NONZERO_PATTERN);
//     totalLeft.AXPY(PETSC_i * laser.getTimeSpacing() / 2.0, Pot, SAME_NONZERO_PATTERN);

//     totalRight.AXPY(-PETSC_i * laser.getTimeSpacing() / 2.0, K, SAME_NONZERO_PATTERN);
//     totalRight.AXPY(-PETSC_i * laser.getTimeSpacing() / 2.0, Pot, SAME_NONZERO_PATTERN);


//     auto firstTermLeft = kroneckerProduct(I,totalLeft,1.0,2*basis.getDegree() + 1);
//     auto firstTermRight = kroneckerProduct(I,totalRight,1.0,2*basis.getDegree() + 1);

//     Matrix modifiedI{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(), 1};
//     for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
//     {   
//         int l{};
//         int m{};
//         std::tie(l,m) = angular.getBlockMap().at(blockIdx);
//         modifiedI.setValue(blockIdx, blockIdx, l*(l+1) / 2.0);
//     }
//     modifiedI.assemble();

//     auto centrifugalLeft = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Invr2,centrifugalLeft,basis,true);
    
//     auto centrifugalRight = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Invr2,centrifugalRight,basis,true);

//     centrifugalLeft *= PETSC_i * laser.getTimeSpacing() / 2.0;
//     centrifugalRight *= -PETSC_i * laser.getTimeSpacing() / 2.0;

//     auto secondTermLeft = kroneckerProduct(modifiedI, centrifugalLeft, 1.0,2*basis.getDegree() + 1);
//     auto secondTermRight = kroneckerProduct(modifiedI, centrifugalRight, 1.0,2*basis.getDegree() + 1);

//     firstTermLeft.AXPY(1.0, secondTermLeft, SAME_NONZERO_PATTERN);
//     firstTermRight.AXPY(1.0, secondTermRight, SAME_NONZERO_PATTERN);

//     return std::make_pair(firstTermLeft,firstTermRight);
// }

// Matrix TDSE::constructZInteraction(const Basis& basis, const Angular& angular)
// {
//     auto Der = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Der,Der,basis,true);


//     auto Invr = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Invr,Invr,basis,true);

//     Matrix Hlm_z_1{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::Z_INT_1,Hlm_z_1,angular);

//     Matrix Hlm_z_2{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::Z_INT_2,Hlm_z_2,angular);
    
//     auto ZInteraction = kroneckerProduct(Hlm_z_1, Der, 2, 2*basis.getDegree() + 1);
//     ZInteraction.AXPY(1.0,kroneckerProduct(Hlm_z_2, Invr, 2, 2*basis.getDegree() + 1),SAME_NONZERO_PATTERN);

//     return ZInteraction;
// }

// std::pair<Matrix,Matrix> TDSE::constructXYInteraction(const Basis& basis, const Angular& angular)
// {
//     auto Der = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Der,Der,basis,true);

//     auto Invr = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Invr,Invr,basis,true);

//     // Compute Hxy_1
//     Matrix Hlm_xy_1{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::XY_INT_1,Hlm_xy_1,angular);

//     Matrix Hlm_xy_2{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::XY_INT_2,Hlm_xy_2,angular);

//     auto Hxy_1 = kroneckerProduct(Hlm_xy_1,Invr,2,2*basis.getDegree() + 1);
//     Hxy_1.AXPY(1.0,kroneckerProduct(Hlm_xy_2,Der,2, 2*basis.getDegree() + 1),SAME_NONZERO_PATTERN);

//     // Construct Hxy_2
//     Matrix Hlm_xy_3{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::XY_INT_3,Hlm_xy_3,angular);

//     Matrix Hlm_xy_4{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),2};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::XY_INT_4,Hlm_xy_4,angular);

//     auto Hxy_2 = kroneckerProduct(Hlm_xy_3,Invr,2,2*basis.getDegree() + 1);
//     Hxy_2.AXPY(1.0,kroneckerProduct(Hlm_xy_4,Der,2,2*basis.getDegree()+1),SAME_NONZERO_PATTERN);

//     return std::make_pair(Hxy_1,Hxy_2);
// }

// Matrix TDSE::constructAtomicS(const Basis& basis, const Angular& angular)
// {
//     Matrix I{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(), 1};
//     for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
//     {
//         I.setValue(blockIdx, blockIdx, 1.0);
//     }
//     I.assemble();

//     auto S = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::S,S,basis,true);
    
//     return kroneckerProduct(I,S,1,2*basis.getDegree() + 1);
// }

// Matrix TDSE::constructXHHG(const Basis& basis, const Angular& angular)
// {
//     auto Invr2 = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Invr2,Invr2,basis,false);

//     Matrix Hlm_hhg_x{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),4};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::X_HHG,Hlm_hhg_x,angular);

//     return kroneckerProduct(Hlm_hhg_x,Invr2,4,2*basis.getDegree() + 1);
// }

// Matrix TDSE::constructYHHG(const Basis& basis, const Angular& angular)
// {
//     auto Invr2 = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Invr2,Invr2,basis,false);

//     Matrix Hlm_hhg_y{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),4};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::Y_HHG,Hlm_hhg_y,angular);

//     return kroneckerProduct(Hlm_hhg_y,Invr2,4,2*basis.getDegree() + 1);
// }

// Matrix TDSE::constructZHHG(const Basis& basis, const Angular& angular)
// {
//     auto Invr2 = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,basis.getNbasis(),basis.getNbasis(),2*basis.getDegree() + 1};
//     RadialMatrix::populateRadialMatrix(RadialMatrixType::Invr2,Invr2,basis,false);

//     Matrix Hlm_hhg_z{PETSC_COMM_SELF,PETSC_DETERMINE,PETSC_DETERMINE,angular.getNlm(),angular.getNlm(),4};
//     AngularMatrix::populateAngularMatrix(AngularMatrixType::Z_HHG,Hlm_hhg_z,angular);

//     return kroneckerProduct(Hlm_hhg_z,Invr2,2,2*basis.getDegree() + 1);
// }

void computeHHG(int rank, std::ofstream& hhgFile, bool hhgStatus, int timeIdx,const Laser& laser, const Vector& state, const Matrix& XHHG, const Matrix& YHHG, const Matrix& ZHHG)
{
    if (hhgStatus)
    {   
        std::complex<double> xVal{};
        std::complex<double> yVal{};
        std::complex<double> zVal{};

        if (laser.getComponents()[0])
        {
            xVal = norm(state,XHHG);
            xVal = xVal * xVal;
        }
        if (laser.getComponents()[1])
        {
            yVal = norm(state,YHHG);
            yVal = yVal  * yVal ;
        }
        if (laser.getComponents()[2])
        {
            zVal = norm(state,ZHHG);
            zVal = zVal * zVal;
        }

        if (rank == 0)
        {
            hhgFile << timeIdx*laser.getTimeSpacing()  << " " << xVal.real() << " "  << laser.A(timeIdx*laser.getTimeSpacing(),0) << " " << yVal.real() << " " << laser.A(timeIdx*laser.getTimeSpacing(),1)   << " " << zVal.real() << " " << laser.A(timeIdx*laser.getTimeSpacing(),2)  << "\n";
        }
    }
}

// void TDSE::solve(int rank,const TISE& tise,const Basis& basis, const Angular& angular, const Atom& atom, const Laser& laser)
// {
//     if (!getStatus())
//     {
//         return;
//     }

//     // Prepate Matrices/Input state for solver
//     auto start_setup = MPI_Wtime();

//     auto initialState = loadInitialState(tise,basis,angular);

//     auto atomicS = constructAtomicS(basis,angular);

//     auto normVal = norm(initialState,atomicS);
//     PetscPrintf(PETSC_COMM_WORLD,"Initial Norm: (%.15f , %.15f) \n",normVal.real(),normVal.imag()); 

//     Matrix interactionLeft{};
//     Matrix interactionRight{};
//     std::tie(interactionLeft,interactionRight) = constructAtomicInteraction(basis,angular,atom,laser);

//     Matrix ZInteraction{};
//     Matrix ZHHG{};
//     if (laser.getComponents()[2])
//     {
//         ZInteraction = constructZInteraction(basis,angular);
//         if (getHHGStatus())
//         {
//             ZHHG = constructZHHG(basis,angular);
//         }
//     }

//     Matrix Hxy_1{};
//     Matrix Hxy_2{};
//     Matrix XHHG{};
//     Matrix YHHG{};
//     if ((laser.getComponents()[0]) || (laser.getComponents()[1]))
//     {
//         std::tie(Hxy_1,Hxy_2) = constructXYInteraction(basis,angular);

//         if (getHHGStatus() && laser.getComponents()[0])
//         {
//             XHHG = constructXHHG(basis,angular);
//         }
//         if (getHHGStatus() && laser.getComponents()[1])
//         {
//             YHHG = constructYHHG(basis,angular);
//         }
//     }

//     PetscScalar alpha = PETSC_i * laser.getTimeSpacing() / 2.0;

//     auto end_setup = MPI_Wtime();
//     PetscPrintf(PETSC_COMM_WORLD,"Time to setup TDSE: %f \n", end_setup-start_setup);

//     // Setup HHG if necessary
//     std::ofstream hhgFile;
//     if (getHHGStatus() && (rank == 0))
//     {
//         hhgFile.open(std::string("misc/") + std::string("hhg_data.txt"));
//         hhgFile << std::fixed << std::setprecision(15); 
//     }


//     // Solve TDSE
//     KSPSolver ksp(PETSC_COMM_WORLD,getMaxIter(),getTol(),getRestart());
//     ksp.setOperators(interactionLeft);

//     auto rhs = Vector{};
//     interactionRight.setupVector(rhs);

//     auto start_solve = MPI_Wtime();
//     for (int timeIdx = 0; timeIdx < laser.getNt(); ++timeIdx)
//     {

//         computeHHG(rank,hhgFile,getHHGStatus(),timeIdx,laser,initialState,XHHG,YHHG,ZHHG);


//         double tNow = timeIdx * laser.getTimeSpacing() + laser.getTimeSpacing() / 2.0;
        
//         if (timeIdx == 0)
//         {
//             if (laser.getComponents()[2])
//             {
//                 auto Az = laser.A(tNow,2);

//                 interactionLeft.AXPY(Az * alpha, ZInteraction,DIFFERENT_NONZERO_PATTERN);
//                 interactionRight.AXPY(-Az * alpha, ZInteraction,DIFFERENT_NONZERO_PATTERN);
//             }
//             if ((laser.getComponents()[0]) || (laser.getComponents()[1]))
//             {
//                 auto deltaAtilde = (laser.A(tNow,0) + PETSC_i * laser.A(tNow,1));
//                 auto deltaAtildeStar = (laser.A(tNow,0) - PETSC_i * laser.A(tNow,1));

//                 interactionLeft.AXPY(alpha * deltaAtildeStar,Hxy_1,DIFFERENT_NONZERO_PATTERN);
//                 interactionRight.AXPY(-alpha * deltaAtildeStar,Hxy_1,DIFFERENT_NONZERO_PATTERN);

//                 interactionLeft.AXPY(alpha*deltaAtilde,Hxy_2,DIFFERENT_NONZERO_PATTERN);
//                 interactionRight.AXPY(-alpha*deltaAtilde,Hxy_2,DIFFERENT_NONZERO_PATTERN);
//             }
//         }
//         else
//         {   
//             double tPrev = (timeIdx - 1) * laser.getTimeSpacing() + laser.getTimeSpacing() / 2.0;

//             if (laser.getComponents()[2])
//             {
//                 auto deltaAz = (laser.A(tNow,2) - laser.A(tPrev,2));

//                 interactionLeft.AXPY( deltaAz * alpha, ZInteraction,SUBSET_NONZERO_PATTERN);
//                 interactionRight.AXPY( -deltaAz * alpha, ZInteraction,SUBSET_NONZERO_PATTERN);
//             }
//             if ((laser.getComponents()[0]) || (laser.getComponents()[1]))
//             {
//                 auto deltaAtilde = ((laser.A(tNow,0) + PETSC_i * laser.A(tNow,1)) - (laser.A(tPrev,0) + PETSC_i * laser.A(tPrev,1)));
//                 auto deltaAtildeStar = ((laser.A(tNow,0) - PETSC_i * laser.A(tNow,1)) - (laser.A(tPrev,0) - PETSC_i * laser.A(tPrev,1)));

//                 interactionLeft.AXPY(alpha * deltaAtildeStar,Hxy_1,SUBSET_NONZERO_PATTERN);
//                 interactionRight.AXPY(-alpha * deltaAtildeStar,Hxy_1,SUBSET_NONZERO_PATTERN);

//                 interactionLeft.AXPY(alpha*deltaAtilde,Hxy_2,SUBSET_NONZERO_PATTERN);
//                 interactionRight.AXPY(-alpha*deltaAtilde,Hxy_2,SUBSET_NONZERO_PATTERN);
//             }
//         }

//         interactionRight.matMult(initialState,rhs);
//         ksp.solve(rhs,initialState);


        
//     }

    
//     auto end_solve = MPI_Wtime();
//     PetscPrintf(PETSC_COMM_WORLD,"Time to solve TDSE: %f \n", end_solve-start_solve);

//     normVal = norm(initialState,atomicS);
//     PetscPrintf(PETSC_COMM_WORLD,"Final Norm: (%.15f , %.15f) \n",normVal.real(),normVal.imag()); 

//     PetscHDF5 viewer(PETSC_COMM_WORLD,getOutputPath(), FILE_MODE_WRITE);
//     viewer.saveVector(outputGroup,outputName,initialState);

//     if (getHHGStatus() && (rank == 0))
//     {
//         hhgFile.close();
//     }

    
// }
