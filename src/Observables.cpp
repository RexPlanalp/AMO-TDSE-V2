#include "Observables.h"

void Block::computeDistribution(int rank,const BSpline& bspline, const TDSE& tdse,const TISE& tise, const Angular& angular)
{
    if (rank != 0)
    {
        return;
    }

    auto S = bspline.PopulateMatrix(PETSC_COMM_SELF,&BSpline::overlapIntegrand, false);

    PetscHDF5 viewer{PETSC_COMM_SELF,tdse.getOutputPath(),FILE_MODE_READ};
    std::string groupName = "";
    std::string vectorName = "psiFinal";
    auto finalState = viewer.loadVector(groupName, vectorName,bspline.getNbasis() * angular.getNlm());

    std::string filename = std::string("misc") + "/block_norms.txt";

    std::ofstream outFile(filename);
    outFile << std::fixed << std::setprecision(15);

    if (!outFile)
    {
        std::cerr << "Error opening file: " << filename << '\n';
    }

    if (getProjOut())
    {
        projectOutBoundStates(finalState,S,tise,angular,bspline);
    }

    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {
        std::cout << "Computing norm for block: " << blockIdx << '\n';

        int start = blockIdx * bspline.getNbasis();

        IndexSet is{PETSC_COMM_SELF,bspline.getNbasis(), start, 1};

        Vector blockVector{};
        VecGetSubVector(finalState.get(), is.get(), &blockVector.get());

        auto normVal = innerProduct(blockVector,S,blockVector);

        VecRestoreSubVector(finalState.get(), is.get(), &blockVector.get());

        outFile << blockIdx << " " << normVal.real() << " " << normVal.imag() << '\n';
    }

    outFile.close();
}

void Block::projectOutBoundStates(Vector& finalState,const Matrix& S,const TISE& tise, const Angular& angular,const BSpline& bspline)
{
    PetscHDF5 viewer(PETSC_COMM_SELF,tise.getOutputPath(),FILE_MODE_READ);

    for (int blockIdx = 0; blockIdx < angular.getNlm(); ++blockIdx)
    {
        int l{};
        int m{};
        std::tie(l,m) = angular.getBlockMap().at(blockIdx);


        int start = blockIdx * bspline.getNbasis();
        
        IndexSet is{PETSC_COMM_SELF,bspline.getNbasis(),start,1};

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
                Vector tiseState = viewer.loadVector(groupName,vectorName,bspline.getNbasis());

                PetscScalar prodVal = innerProduct(tiseState,S,blockVector);

                blockVector.AXPY(-prodVal, tiseState);
            }


        }

        VecRestoreSubVector(finalState.get(),is.get(),&blockVector.get());
    }
}