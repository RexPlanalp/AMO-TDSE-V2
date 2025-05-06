#include "Box.h"

int Box::getNr() const 
{
    return static_cast<int>(std::floor(getGridSize() / getGridSpacing()));
}

double Box::getPosition(int i) const
{
    return i * getGridSpacing();
}