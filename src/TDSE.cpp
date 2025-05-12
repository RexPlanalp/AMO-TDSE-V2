#include "TDSE.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscMat.h"
#include "PetscWrappers/PetscOperators.h"

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
        std::cout << "Krylov Dimensionality: " << getKrylovDim() << "\n\n";
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

void TDSE::solve(const TISE& tise,const BSpline& bspline, const Angular& angular, const Atom& atom, const Laser& laser)
{
    auto initialState = loadInitialState(tise,bspline,angular);

    Matrix interactionLeft{};
    Matrix interactionRight{};

    std::tie(interactionLeft,interactionRight) = constructAtomicInteraction(bspline,angular,atom,laser);

    
}