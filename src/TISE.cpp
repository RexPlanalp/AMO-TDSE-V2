#include "TISE.h"

void TISE::validateInput()
{
    if (MAXITER() <= 0)
    {
        throw std::runtime_error("max_iter must be greater than zero. You entered: " + std::to_string(MAXITER()));
    }
    if (NMAX() <= 0)
    {
        throw std::runtime_error("nmax must be greather than zero. You entered: " + std::to_string(NMAX()));
    }
}