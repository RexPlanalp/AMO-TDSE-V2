#pragma once
#include "petscmat.h"

class Matrix 
{
    public:
        
        // Default Constructor
        Matrix() : m{nullptr}, iStart{0}, iEnd{0} {}

        // Regular Constructor
        Matrix(MPI_Comm comm, PetscInt localRows, PetscInt localCols, PetscInt globalRows, PetscInt globalCols, PetscInt nnz_p)
        {
            MatCreate(comm, &m);
            MatSetSizes(m , localRows, localCols, globalRows, globalCols);
            if (comm == PETSC_COMM_WORLD)
            {
                MatMPIAIJSetPreallocation(m, nnz_p, NULL, nnz_p, NULL);
            }
            else if (comm == PETSC_COMM_SELF)
            {
                MatSeqAIJSetPreallocation(m, nnz_p, NULL);
            }
            MatSetFromOptions(m);
            MatGetOwnershipRange(m, &iStart, &iEnd);
        }
        
        // Destructor 
        ~Matrix()
        {
            MatDestroy(&m);
        }

        // Copy Constructor
        Matrix(const Matrix& other)
        {
            MatDuplicate(other.m, MAT_COPY_VALUES, &m);
        }

        // Copy Assignment
        Matrix& operator=(Matrix other) noexcept
        {
            swap(*this, other);
            return *this;
        }

        // Move Constructor
        Matrix(Matrix&& other) noexcept
        : m{other.m}
        {
            other.m = nullptr;
            other.iStart = 0;
            other.iEnd = 0;
        }
        
        // Move Assignment 
        Matrix& operator=(Matrix&& other) noexcept
        {
            if (this != &other)
            {
                if (m)
                {
                    MatDestroy(&m);
                }
                m = other.m;
                other.m = nullptr;
                other.iStart = 0;
                other.iEnd = 0;
            }
            return *this;
        }
        friend void swap(Matrix& a, Matrix& b) noexcept
        {
            std::swap(a.m,b.m);
            std::swap(a.iStart, b.iStart);
            std::swap(a.iEnd,   b.iEnd);
        }
    
    void assemble()
    {
        MatAssemblyBegin(get(), MAT_FINAL_ASSEMBLY); 
        MatAssemblyEnd(get(), MAT_FINAL_ASSEMBLY); 
    }

    void AXPY(PetscScalar alpha, const Matrix& other, MatStructure pattern)
    {
        MatAXPY(get(), alpha,other.get(),pattern);
    }

    Mat& get() { return m; }
    const Mat& get() const {return m;}

    PetscInt getStart() { return iStart;}
    PetscInt getEnd() { return iEnd;}


    private:
        Mat m;
        PetscInt iStart{};
        PetscInt iEnd{};
};