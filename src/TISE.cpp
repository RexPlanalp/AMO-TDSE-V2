#include "TISE.h"
#include "Basis.h"
#include "PetscWrappers/PetscEPS.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscOperators.h"
#include "Angular.h"
#include "Atom.h"
#include "MatrixElements.h"
#include "common.h"




// void TISE::printConfiguration(int rank)
// {
//     if (rank == 0)
//     {
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
//         std::cout << "TISE Configuration: " << "\n\n";
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
//         std::cout << "status: " << getStatus() << "\n\n";
//         std::cout << "maxIter: " << getMaxIter() << "\n\n";
//         std::cout << "outputPath: " << getOutputPath() << "\n\n";
//         std::cout << "nmax: " << getNmax() << "\n\n";
//     }
// }
