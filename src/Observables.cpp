#include "Observables.h"

#include <cmath>


void Observables::buildNe()
{
    m_Ne = static_cast<int>(std::round(getEmax() / getEmin())) + 1;
}

