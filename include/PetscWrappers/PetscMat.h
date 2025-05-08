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
            if (comm = PETSC_COMM_WORLD)
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
        Matrix& operator=(Matrix other)
        {
            swap(*this, other);
            return *this;
        }

        // Move Constructor
        Matrix(Matrix&& other) noexcept
        : m{other.m}
        {
            other.m = nullptr;
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
            }
            return *this;
        }
        friend void swap(Matrix& a, Matrix& b) noexcept
        {
        Mat tmp = a.m;
        a.m = b.m;
        b.m = tmp;
        }
    

    Mat& get() { return m; }
    const Mat& get() const {return m;}


    private:
        Mat m;
        PetscInt iStart{};
        PetscInt iEnd{};
};