#include "Atom.h"

void Atom::validateInput()
{
    if (species == std::string{})
    {
        throw std::invalid_argument("Species must be non-empty. You entered: " + species);
    }
}