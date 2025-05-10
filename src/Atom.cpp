#include "Atom.h"

void Atom::printConfiguration(int rank)
{
    if (rank == 0)
    {
        std::cout << "Atom Configuration: " << "\n\n";
        std::cout << "Potential: " << species <<  "\n\n";
    }
}