#include "Observables.h"

void Block::computeDistribution(int rank,const BSpline& bspline, const TDSE& tdse, const Angular& angular)
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