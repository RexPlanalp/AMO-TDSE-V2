#pragma once
#include "petscvec.h"




class Vector
{
  public:

    // Default Constructor
    Vector() : v{nullptr}, iStart{0}, iEnd{0} {}
    
    // Regular Constructor
    Vector(MPI_Comm comm, PetscInt localSize, PetscInt globalSize)
    {
      VecCreate(comm, &v);
      VecSetSizes(v,localSize,globalSize);
      if (comm == PETSC_COMM_WORLD)
      {
        VecSetType(v, VECMPI);
      }
      else if (comm == PETSC_COMM_SELF)
      {
        VecSetType(v, VECSEQ);
      }
      VecSetFromOptions(v);
      VecGetOwnershipRange(v,&iStart,&iEnd);
    }

    // Destructor
    ~Vector()
    {
      if (v)
      {
        VecDestroy(&v);
      }
    }

    // Copy Constructor
    Vector(const Vector& other)
    {
      VecDuplicate(other.v, &v);
      VecCopy(other.v, v);
    }

    // Copy Assignment
    Vector& operator=(Vector other) noexcept
    {
      swap(*this, other);
      return *this;
    }

    // Move Constructor
    Vector(Vector&& other) noexcept
    : v{other.v}
    {
      other.v = nullptr;
      other.iStart = 0;
      other.iEnd= 0;
    }

    // Move Assignment
    Vector& operator=(Vector&& other) noexcept
    {
      if (this != &other)
      {
        if (v) 
        {
          VecDestroy(&v);
        }
        v = other.v;
        other.v = nullptr;
        other.iStart = 0;
        other.iEnd = 0;
      }
      return *this;
    }

    friend void swap(Vector& a, Vector& b) noexcept
    {
      std::swap(a.v,b.v);
      std::swap(a.iStart, b.iStart);
      std::swap(a.iEnd,   b.iEnd);
    }

    Vec& get() { return v; }
    const Vec& get() const {return v;}

    PetscInt getStart() const { return iStart;}
    PetscInt getEnd() const { return iEnd;}

    

    void setValue(PetscInt idx, PetscScalar value)
    {
      VecSetValue(get(), idx, value, INSERT_VALUES);
    }

    void assemble()
    {
      VecAssemblyBegin(get());
      VecAssemblyEnd(get());
    }


    Vector& operator*=(PetscScalar alpha) 
    {
      VecScale(v, alpha);
      return *this;
    }


  private:
    Vec v;
    PetscInt iStart{};
    PetscInt iEnd{};
};
