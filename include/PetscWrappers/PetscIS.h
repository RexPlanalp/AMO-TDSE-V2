#pragma once
#include "petscis.h"




class IndexSet
{
  public:

    // Default Constructor
    IndexSet() : is{nullptr} {}
    
    // Regular Constructor
    IndexSet(MPI_Comm comm, PetscInt size, PetscInt start, PetscInt stepSize)
    {
      ISCreateStride(comm,size,start,stepSize,&is);
    }

    // Destructor
    ~IndexSet()
    {
      if (is)
      {
        ISDestroy(&is);
      }
    }

    const IS& get() const {return is;}
    IS& get() {return is;}


  private:
    IS is;
};
