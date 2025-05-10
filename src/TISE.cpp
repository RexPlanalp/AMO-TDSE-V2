#include "TISE.h"
#include "BSpline.h"
#include "PetscWrappers/PetscEPS.h"
#include "Angular.h"

void TISE::validateInput()
{
    if (MAXITER() <= 0)
    {
        throw std::runtime_error("max_iter must be greater than zero. You entered: " + std::to_string(MAXITER()));
    }
    if (NMAX() <= 0)
    {
        throw std::runtime_error("nmax must be greather than zero. You entered: " + std::to_string(NMAX()));
    }
}

void TISE::solve(const BSpline& bspline, const Atom& atom, const Angular& angular)
{
    Matrix K = bspline.PopulateMatrix([&](int ii, int jj, std::complex<double> x)
    {
        return bspline.kineticIntegrand(ii, jj, x);
    }
    , false);

    Matrix Invr2 = bspline.PopulateMatrix([&](int ii, int jj, std::complex<double> x)
    {
        return bspline.invr2Integrand(ii, jj, x);
    }
    , false);

    
    Matrix Pot = bspline.PopulateMatrix([&](int ii, int jj, std::complex<double> x)
    {
        return bspline.HIntegrand(ii, jj, x);
    }
    , false); 
    

    Matrix S = bspline.PopulateMatrix([&](int ii, int jj, std::complex<double> x)
    {
        return bspline.overlapIntegrand(ii, jj, x);
    }
    , false);

    K.AXPY(1.0, Pot, SAME_NONZERO_PATTERN);
    

    EPSSolver epssolver{PETSC_COMM_WORLD,MAXITER(),Tol()};
   

    for (int l{0}; l < angular.LMax(); ++l)
    {
        if (l > 0)
        {
            K.AXPY(l, Invr2,SAME_NONZERO_PATTERN);
        }

        int pairs = NMAX() - l;
        epssolver.reset();
        epssolver.setOperators(K,S);
        epssolver.setDimensions(pairs);
        epssolver.solve();

        for (int idx{0}; idx < epssolver.NCONV(); ++idx)
        {
            PetscScalar eigenvalue = epssolver.getEigenvalue(idx);
            PetscPrintf(
                PETSC_COMM_WORLD,
                "Eigenvalue [%d] = %g + %g i\n",
                (int)idx,
                PetscRealPart(eigenvalue),
                PetscImaginaryPart(eigenvalue)
              );
        }
    }






   
}