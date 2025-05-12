#include "TDSE.h"


Matrix KroneckerProduct(const Matrix& A, const Matrix& B,PetscInt nnz_A, PetscInt nnz_B) 
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
        std::cout << "Initial State: " << "n = " << getInitialState()[0]  << " l = " << getInitialState()[1]  << " m = " << getInitialState()[2] << "\n\n";

    }
}