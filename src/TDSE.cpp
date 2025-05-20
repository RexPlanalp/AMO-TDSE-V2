#include "TDSE.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscMat.h"
#include "PetscWrappers/PetscOperators.h"
#include "PetscWrappers/PetscKSP.h"
#include "MatrixElements.h"

#include "MatrixElements.h"


// void TDSE::printConfiguration(int rank)
// {
//     if (rank == 0)
//     {
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
//         std::cout << "TDSE Configuration: " << "\n\n";
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
//         std::cout << "status: " << getStatus() << "\n\n";
//         std::cout << "maxIter: " << getMaxIter() << "\n\n";
//         std::cout << "outputPath: " << getOutputPath() << "\n\n";
//         std::cout << "Restart: " << getRestart() << "\n\n";
//         std::cout << "Initial nlm: " << "n = " << getInitialNLM()[0]  << " l = " << getInitialNLM()[1]  << " m = " << getInitialNLM()[2] << "\n\n";

//     }
//}














