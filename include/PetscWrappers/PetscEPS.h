#pragma once
#include "slepceps.h"

#include "PetscWrappers/PetscMat.h"
#include "PetscWrappers/PetscVec.h"

class EPSSolver
{
  public:

    EPSSolver(MPI_Comm comm, PetscInt maxIter, PetscReal tolerance)
    : maxIter{maxIter}
    , tolerance{tolerance}
    {
      EPSCreate(comm, &eps);
      EPSSetProblemType(eps,EPS_GNHEP);
      EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);
      EPSSetType(eps,EPSKRYLOVSCHUR);
      EPSSetTolerances(eps,tolerance,maxIter);
      EPSSetFromOptions(eps);
    }



    ~EPSSolver()
    {
      if (eps) 
      {
        EPSDestroy(&eps);
      }
    }


    PetscInt getNconv() const {return nconv;}

    void setOperators(const Matrix& H, const Matrix& S)
    {
        EPSSetOperators(eps,H.get(),S.get());
    }

    void setDimensions(PetscInt pairs)
    {
      EPSSetDimensions(eps,pairs, PETSC_DEFAULT,PETSC_DEFAULT);
    }

    void solve()
    {
      EPSSolve(eps);
      EPSGetConverged(eps,&nconv);
    }

    PetscScalar getEigenvalue(PetscInt i)
    {
      PetscScalar eigenvalue{};
      EPSGetEigenvalue(eps, i, &eigenvalue, nullptr);
      return eigenvalue;
    }


    void getEigenvector(PetscInt i, Vector& eigenvector)
    {
      EPSGetEigenvector(eps, i ,eigenvector.get(),nullptr);
    }

    const EPS& get() const {return eps;}
    EPS& get() {return eps;}

    


  private:
    EPS eps{nullptr};

    // Member List Initialized
    PetscInt maxIter{};
    PetscReal tolerance{};

    // Default initialized
    PetscInt nconv;
    
};
