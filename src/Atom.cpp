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

