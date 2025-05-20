#include "Atom.h"


void Atom::buildPotential()
{
    if (getSpecies()== "H") 
    {
        m_potential  = Potentials::hydrogen;
        m_derivative = Potentials::hydrogenDeriv;
        m_type = RadialMatrixType::H;
    }
}

// void Atom::printConfiguration(int rank)
// {
//     if (rank == 0)
//     {   
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
//         std::cout << "Atom Configuration: " << "\n\n";
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";

//         std::cout << "Species: " << species <<  "\n\n";
//     }
// }