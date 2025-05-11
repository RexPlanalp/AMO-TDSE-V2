#pragma once
#include "slepceps.h"

#include "PetscWrappers/PetscMat.h"
#include "PetscWrappers/PetscVec.h"

class EPSSolver
{
  public:

    EPSSolver(MPI_Comm comm, PetscInt maxIter, PetscReal tolerance)
    : maxIter{maxIter}, tolerance{tolerance}
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

    EPSSolver(const EPSSolver&)            = delete;
    EPSSolver& operator=(const EPSSolver&) = delete;
    EPSSolver(EPSSolver&&)                 = delete;
    EPSSolver& operator=(EPSSolver&&)      = delete;

    PetscInt getNconv() const {return nconv;}


    void reset()
    {
        EPSReset(eps);
    }

    void setOperators(const Matrix& H, const Matrix& S)
    {
        EPSSetOperators(eps,H.get(),S.get());
    }

    void setDimensions(PetscInt pairs)
    {
      EPSSetDimensions(eps,pairs, PETSC_DETERMINE ,PETSC_DETERMINE);
    }

    void solve()
    {
      EPSSolve(eps);
      EPSGetConverged(eps,&nconv);
    }

    PetscScalar getEigenvalue(PetscInt i)
    {
      PetscScalar eigenvalue{};
      EPSGetEigenvalue(eps, i, &eigenvalue, NULL);
      return eigenvalue;
    }


    Vector getEigenvector(PetscInt i, const Matrix& S)
    {
      Vector eigenvector{};
      MatCreateVecs(S.get(), &eigenvector.get(),NULL);
      EPSGetEigenvector(eps, i ,eigenvector.get(),NULL);

      return eigenvector;

    }

    


  private:
    EPS eps{nullptr};

    // Member List Initialized
    PetscInt maxIter{};
    PetscReal tolerance{};

    // Default initialized
    PetscInt nconv;
    
};
