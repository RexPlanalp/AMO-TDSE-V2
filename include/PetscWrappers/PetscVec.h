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
    Vector& operator=(Vector other)
    {
      swap(*this, other);
      return *this;
    }

    // Move Constructor
    Vector(Vector&& other) noexcept
    : v{other.v}
    {
      other.v = nullptr;
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
      }
      return *this;
    }

    friend void swap(Vector& a, Vector& b) noexcept
    {
      Vec tmp = a.v;
      a.v = b.v;
      b.v = tmp;
    }

    Vec& get() { return v; }
    const Vec& get() const {return v;}

    PetscInt getStart() { return iStart;}
    PetscInt getEnd() { return iEnd;}


  private:
    Vec v;
    PetscInt iStart{};
    PetscInt iEnd{};
};
