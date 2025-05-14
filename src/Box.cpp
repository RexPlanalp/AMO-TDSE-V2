#include "Box.h"

#include "common.h"

double Box::getPosition(int i) const
{
    return positions[i];
}

void Box::printConfiguration(int rank) const
{
    if (rank == 0)
    {
        std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        std::cout << "Box Configuration: " << "\n\n";
        std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
        std::cout << "rmax: " << getGridSize() << "\n\n";
        std::cout << "dr: " << getGridSpacing() << "\n\n";
        std::cout << "Nr: " << getNr() << "\n\n";
    }
}