#pragma once

#include "common.h"

namespace BSplines
{
    inline std::complex<double> dB(int i,int degree, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) 
    {
        const int p = degree;
        if (p == 0)
        {
            return std::complex<double>(0.0, 0.0);
        }
        
        auto compute_basis = [&](int ii) 
        {
            std::vector<std::complex<double>> N(p);

            for (int j = 0; j < p; ++j) 
            {
                double t0 = std::real(localKnots[ii + j]);
                double t1 = std::real(localKnots[ii + j + 1]);
                N[j] = (t0 <= x.real() && x.real() < t1)
                    ? std::complex<double>(1.0, 0.0)
                    : std::complex<double>(0.0, 0.0);
            }

            for (int d = 1; d < p; ++d) 
            {
                for (int j = 0; j < p - d; ++j) 
                {
                    auto ld = localKnots[ii + j + d]     - localKnots[ii + j];
                    auto rd = localKnots[ii + j + d + 1] - localKnots[ii + j + 1];
                    std::complex<double> a(0.0, 0.0), b(0.0, 0.0);
                    if (std::abs(ld) > 0.0)
                        a = (x - localKnots[ii + j])       / ld * N[j];
                    if (std::abs(rd) > 0.0)
                        b = (localKnots[ii + j + d + 1] - x) / rd * N[j + 1];
                    N[j] = a + b;
                }
            }
            return N[0];
        };

        std::complex<double> N_i   = compute_basis(i);
        std::complex<double> N_ip1 = compute_basis(i + 1);

        auto denom1 = localKnots[i + p]     - localKnots[i];
        auto denom2 = localKnots[i + p + 1] - localKnots[i + 1];
        double pd = static_cast<double>(p);

        std::complex<double> term1(0.0, 0.0), term2(0.0, 0.0);
        if (std::abs(denom1) > 0.0)
            term1 = pd / denom1 * N_i;
        if (std::abs(denom2) > 0.0)
            term2 = pd / denom2 * N_ip1;

        return term1 - term2;
    }

    inline std::complex<double> B(int i, int degree, std::complex<double> x, const std::vector<std::complex<double>>& localKnots)
    {
        std::vector<std::complex<double>> N(degree + 1);

        for (int j = 0; j <= degree; ++j) 
        {
            double t0 = std::real(localKnots[i + j]);
            double t1 = std::real(localKnots[i + j + 1]);
            N[j] = (t0 <= std::real(x) && std::real(x) < t1)
                ? std::complex<double>(1.0, 0.0)
                : std::complex<double>(0.0, 0.0);
        }

        for (int d = 1; d <= degree; ++d) 
        {
            for (int j = 0; j <= degree - d; ++j) 
            {
                auto left_den  = localKnots[i + j + d]     - localKnots[i + j];
                auto right_den = localKnots[i + j + d + 1] - localKnots[i + j + 1];

                std::complex<double> t1(0.0, 0.0), t2(0.0, 0.0);

                if (std::abs(left_den)  > 0.0)
                    t1 = (x - localKnots[i + j])       / left_den  * N[j];
                if (std::abs(right_den) > 0.0)
                    t2 = (localKnots[i + j + d + 1] - x) / right_den * N[j + 1];

                N[j] = t1 + t2;
            }
        }

        return N[0];
    }

}