#include "Box.h"

void Box::validateInput()
{
    if (grid_size <= 0.0)
    {
        throw std::invalid_argument("Grid size must be greater than zero. You entered: " + std::to_string(GridSize()));
    }
    if (grid_spacing <= 0.0)
    {
        throw std::invalid_argument("Grid spacing must be greater than zero. You entered: " + std::to_string(GridSpacing()));
    }
}

int Box::Nr() const 
{
    return static_cast<int>(std::round(GridSize() / GridSpacing())) + 1;
}

double Box::Position(int i) const
{
    return i * GridSpacing();
}