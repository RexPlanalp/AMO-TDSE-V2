#include "TISE.h"
#include "BSpline.h"
#include "PetscWrappers/PetscEPS.h"
#include "Angular.h"
#include "Atom.h"

#include "common.h"

void TISE::solve(const BSpline& bspline, const Atom& atom, const Angular& angular)
{   
    if (!getStatus())
    {
        return;
    }
   

    Matrix K = bspline.PopulateMatrix(&BSpline::kineticIntegrand,false);

    Matrix Invr2 = bspline.PopulateMatrix(&BSpline::invr2Integrand,false);
    
    Matrix Pot = bspline.PopulateMatrix(&BSpline::HIntegrand,false);
    

    Matrix S = bspline.PopulateMatrix(&BSpline::overlapIntegrand,false);

    K.AXPY(1.0, Pot, SAME_NONZERO_PATTERN);

    
    

    EPSSolver epssolver{PETSC_COMM_WORLD,getMaxIter(),getTol()};
   

    for (int l{0}; l < angular.getLmax(); ++l)
    {
        if (l > 0)
        {
            K.AXPY(l, Invr2,SAME_NONZERO_PATTERN);
        }

        int pairs = getNmax() - l;
        if (pairs <=0)
        {
            continue;
        }   

        

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