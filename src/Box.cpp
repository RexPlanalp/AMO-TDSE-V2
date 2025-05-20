#include "Box.h"

#include "common.h"

double Box::operator[](int i) const
{
    return m_positions[i];
}

void Box::buildNr()
{
     m_Nr = static_cast<int>(std::round(getGridSize() / getGridSpacing())) + 1;
}

void  Box::buildPositions()
{
    m_positions.resize(getNr());
    for (int idx = 0; idx < getNr(); ++idx)
    {
        m_positions[idx] = idx * getGridSpacing();
    }
}


