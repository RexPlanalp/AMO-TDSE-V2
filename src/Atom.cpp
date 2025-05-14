#include "Atom.h"

void Atom::printConfiguration(int rank)
{
    if (rank == 0)
    {   
        std::cout << std::setfill('\\') << std::setw(24) << "" << '\n';
        std::cout << "Atom Configuration: " << "\n\n";
        std::cout << std::setfill('\\') << std::setw(24) << "" << '\n';

        std::cout << "Species: " << species <<  "\n\n";
    }
}