#pragma once

#include <complex>

namespace Potentials
{
    inline std::complex<double> hydrogen(const std::complex<double>& x) 
    {
        return -std::complex<double>{1.0} / (x + std::complex<double>{1e-25});
    }

    inline std::complex<double> hydrogenDeriv(const std::complex<double>& x) 
    {
        return  std::complex<double>{1.0} / (x * x + std::complex<double>{1e-25});
    }
}