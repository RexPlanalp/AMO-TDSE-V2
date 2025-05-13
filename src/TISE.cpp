#include "TISE.h"
#include "BSpline.h"
#include "PetscWrappers/PetscEPS.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscOperators.h"
#include "Angular.h"
#include "Atom.h"

#include "common.h"

void TISE::solve(const BSpline& bspline, const Atom& atom, const Angular& angular)
{   
    if (!getStatus())
    {
        return;
    }
   

    Matrix K = bspline.PopulateMatrix(PETSC_COMM_WORLD,&BSpline::kineticIntegrand,false);

    Matrix Invr2 = bspline.PopulateMatrix(PETSC_COMM_WORLD,&BSpline::invr2Integrand,false);
    
    Matrix Pot{};

    if (atom.getSpecies() == "H")
    {
        Pot = bspline.PopulateMatrix(PETSC_COMM_WORLD,&BSpline::HIntegrand,false);
    }
     
    

    Matrix S = bspline.PopulateMatrix(PETSC_COMM_WORLD,&BSpline::overlapIntegrand,false);

    K.AXPY(1.0, Pot, SAME_NONZERO_PATTERN);

    
    
    PetscHDF5 viewer{PETSC_COMM_WORLD, getOutputPath(), FILE_MODE_WRITE};
    EPSSolver epssolver{PETSC_COMM_WORLD,getMaxIter(),getTol()};
    epssolver.setOperators(K,S);
    std::string eigenvalueGroup = "eigenvalues";
    std::string eigenvectorGroup = "eigenvectors";
   

    for (int l = 0; l < angular.getLmax(); ++l)
    {
        if (l > 0)
        {
            K.AXPY(l, Invr2,SAME_NONZERO_PATTERN);
        }

        int reqPairs = getNmax() - l;
        if (reqPairs <=0)
        {
            continue;
        }   

        PetscPrintf(PETSC_COMM_WORLD,"Solving for l = %d \n\n",l); 
        //epssolver.setOperators(K,S);
        epssolver.setDimensions(reqPairs);

        auto start = MPI_Wtime();
        epssolver.solve();
        auto end = MPI_Wtime();
        PetscPrintf(PETSC_COMM_WORLD,"Solved eigenvalue problem in:  %f seconds \n", end-  start);

      
        for (int convPairs = 0; convPairs < epssolver.getNconv(); ++convPairs)
        {
            auto eigenvalue = epssolver.getEigenvalue(convPairs);
            auto eigenvector = epssolver.getEigenvector(convPairs,S);

            normalize(eigenvector,S);

            if (PetscRealPart(eigenvalue) > 0)
            {
                continue;
            }

            std::string eigenvectorName = std::string("psi_") + std::to_string(convPairs + l + 1) + std::string("_") + std::to_string(l);
            std::string eigenvalueName = std::string("E_") + std::to_string(convPairs + l + 1) + std::string("_") + std::to_string(l);

            viewer.saveValue(eigenvalueGroup, eigenvalueName, eigenvalue);
            viewer.saveVector(eigenvectorGroup,eigenvectorName, eigenvector);


            auto normVal = norm(eigenvector,S);
            PetscPrintf(PETSC_COMM_WORLD,"Eigenvector %d -> Norm(%.4f , %.4f) -> Eigenvalue(%.4f , %.4f)  \n",convPairs+1,normVal.real(),normVal.imag(),eigenvalue.real(),eigenvalue.imag()); 
        }
        PetscPrintf(PETSC_COMM_WORLD,"\n");
        
    }

}


void TISE::printConfiguration(int rank)
{
    if (rank == 0)
    {
        std::cout << "TISE Configuration: " << "\n\n";
        std::cout << "status: " << getStatus() << "\n\n";
        std::cout << "maxIter: " << getMaxIter() << "\n\n";
        std::cout << "outputPath: " << getOutputPath() << "\n\n";
        std::cout << "nmax: " << getNmax() << "\n\n";
    }
}
