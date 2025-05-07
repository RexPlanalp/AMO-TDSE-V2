#include "Box.h"

int Box::Nr() const 
{
    return static_cast<int>(std::floor(GridSize() / GridSpacing())) + 1;
}

double Box::Position(int i) const
{
    return i * GridSpacing();
}