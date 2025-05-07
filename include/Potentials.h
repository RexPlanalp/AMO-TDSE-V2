#pragma once


namespace Potentials
{
    template<typename T>
    inline T hydrogenPotential(T position)
    {
        return - 1.0 / (position + 1E-25);
    }

    template<typename T>
    inline T hydrogenPotentialDerivative(T position)
    {
        return 1.0 / (position * position + 1E-25);
    }
}
