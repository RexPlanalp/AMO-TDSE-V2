#include "Observables.h"

void Observable::computeDistribution(int rank,const Basis& Basis, const TDSE& tdse,const TISE& tise, const Angular& angular)
{
    if (rank != 0)
    {
        return;
    }

    auto S = Matrix{PETSC_COMM_SELF,PETSC_DECIDE,PETSC_DECIDE,Basis.getNbasis(),Basis.getNbasis(),2*Basis.getDegree() + 1};
    RadialMatrix::populateRadialMatrix(RadialMatrixType::S,S,Basis,false);

    PetscHDF5 viewer{PETSC_COMM_SELF,tdse.getOutputPath(),FILE_MODE_READ};
    std::string groupName = "";
    std::string vectorName = "psiFinal";
    auto finalState = viewer.loadVector(groupName, vectorName,Basis.getNbasis() * angular.getNlm());

    std::string filename = std::string("misc") + "/block_norms.txt";

    std::ofstream outFile(filename);
    outFile << std::fixed << std::setprecision(15);

    if (!outFile)
    {
        std::cerr << "Error opening file: " << filename << '\n';
    }

    if (getProjOut())
    {
        projectOutBoundStates(finalState,S,tise,angular,Basis);
    }

    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {
        std::cout << "Computing norm for block: " << blockIdx << '\n';

        int start = blockIdx * Basis.getNbasis();

        IndexSet is{PETSC_COMM_SELF,Basis.getNbasis(), start, 1};

        Vector blockVector{};
        VecGetSubVector(finalState.get(), is.get(), &blockVector.get());

        auto normVal = innerProduct(blockVector,S,blockVector);

        VecRestoreSubVector(finalState.get(), is.get(), &blockVector.get());

        outFile << blockIdx << " " << normVal.real() << " " << normVal.imag() << '\n';
    }

    outFile.close();
}

void Observable::projectOutBoundStates(Vector& finalState,const Matrix& S,const TISE& tise, const Angular& angular,const Basis& Basis)
{
    PetscHDF5 viewer(PETSC_COMM_SELF,tise.getOutputPath(),FILE_MODE_READ);

    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockIdx);


        int start = blockIdx * Basis.getNbasis();
        
        IndexSet is{PETSC_COMM_SELF,Basis.getNbasis(),start,1};

        auto blockVector = Vector{};

        VecGetSubVector(finalState.get(),is.get(),&blockVector.get());


        for (int n = 0; n <= tise.getNmax(); ++n)
        {
            std::string groupName = "eigenvectors";
            std::string vectorName = std::string("psi_") + std::to_string(n) + std::string("_") + std::to_string(l);

            std::string datasetName = groupName + "/" + vectorName;

            PetscBool hasDataset{};
            PetscViewerHDF5HasDataset(viewer.get(),datasetName.c_str(),&hasDataset);

            if (hasDataset)
            {
                Vector tiseState = viewer.loadVector(groupName,vectorName,Basis.getNbasis());

                PetscScalar prodVal = innerProduct(tiseState,S,blockVector);

                blockVector.AXPY(-prodVal, tiseState);
            }


        }

        VecRestoreSubVector(finalState.get(),is.get(),&blockVector.get());
    }
}