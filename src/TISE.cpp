#include "TISE.h"
#include "Basis.h"
#include "PetscWrappers/PetscEPS.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscOperators.h"
#include "Angular.h"
#include "Atom.h"

#include "common.h"

void TISE::solve(const Basis& Basis, const Atom& atom, const Angular& angular)
{   
    // If we are not running the TISE, exit
    if (!getStatus())
    {
        return;
    }

    // Start measuring at the very beginning
    auto total_start = MPI_Wtime();

    // Create HDF5 File to store output
    PetscHDF5 viewer{PETSC_COMM_WORLD, getOutputPath(), FILE_MODE_WRITE};
    
    // Create Kinetic Matrix
    Matrix K = Basis.PopulateMatrix(PETSC_COMM_WORLD,&Basis::kineticIntegrand,false);

    // Create Inverse r squared matrix
    Matrix Invr2 = Basis.PopulateMatrix(PETSC_COMM_WORLD,&Basis::invr2Integrand,false);
    
    // Create 
    Matrix Pot{};

    if (atom.getSpecies() == "H")
    {
        Pot = Basis.PopulateMatrix(PETSC_COMM_WORLD,&Basis::HIntegrand,false);
    }
     
    Matrix S = Basis.PopulateMatrix(PETSC_COMM_WORLD,&Basis::overlapIntegrand,false);

    K.AXPY(1.0, Pot, SAME_NONZERO_PATTERN);

    EPSSolver epssolver{PETSC_COMM_WORLD,getMaxIter(),getTol()};
    epssolver.setOperators(K,S);

    auto eigenvector = Vector{};
    S.setupVector(eigenvector);
    
    for (int l = 0; (l < angular.getLmax()) && (l < getNmax()); ++l)
    {
        if (l > 0)
        {
            K.AXPY(l, Invr2,SAME_NONZERO_PATTERN);
        }
        
        PetscPrintf(PETSC_COMM_WORLD,"Solving for l = %d \n",l); 

        int reqPairs = getNmax() - l;
        epssolver.setDimensions(reqPairs);
        epssolver.solve();

        for (int convPair = 0; convPair < epssolver.getNconv(); ++convPair)
        {
            auto eigenvalue = epssolver.getEigenvalue(convPair);
            epssolver.getEigenvector(convPair,eigenvector);

            normalize(eigenvector,S);

            if (PetscRealPart(eigenvalue) > 0)
            {
                continue;
            }

            std::string eigenvectorName = std::string("psi_") + std::to_string(convPair + l + 1) + std::string("_") + std::to_string(l);
            std::string eigenvalueName = std::string("E_") + std::to_string(convPair + l + 1) + std::string("_") + std::to_string(l);

            viewer.saveValue(eigenvalueGroup, eigenvalueName, eigenvalue);
            viewer.saveVector(eigenvectorGroup,eigenvectorName, eigenvector);

            auto normVal = norm(eigenvector,S);
            PetscPrintf(PETSC_COMM_WORLD,"Eigenvector %d -> Norm(%.4f , %.4f) -> Eigenvalue(%.4f , %.4f)  \n",convPair+1,normVal.real(),normVal.imag(),eigenvalue.real(),eigenvalue.imag()); 
        }
        auto total_end = MPI_Wtime();
        PetscPrintf(PETSC_COMM_WORLD,"Solved eigenvalue problem in:  %f seconds \n\n", total_end -  total_start);
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
