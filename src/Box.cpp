#include "Box.h"

int Box::Nr() const 
{
    return static_cast<int>(std::round(GridSize() / GridSpacing())) + 1;
}

double Box::Position(int i) const
{
    return i * GridSpacing();
}