
#include "BSpline.h"
#include <complex>
#include <fstream>
#include <iostream>
#include "Box.h"
#include "GaussLegendre.h"

void BSpline::validateInput()
{
    if (NBSpline() <= 0)
    {
        throw std::invalid_argument("n_bspline must be greater than zero. You entered: " + std::to_string(NBSpline()));
    }
    if (Order() < 0)
    {
        throw std::invalid_argument("order must be greater than or equal to zero. You entered: " + std::to_string(Order()));
    }
    if (!((R0_R() >= 0.0) && (R0_R() <= 1.0)))
    {
        throw std::invalid_argument("R0_r must be between 0.0 and 1.0. You entered: " + std::to_string(R0_R()));
    }
    if (!((ETA_R() >= 0.0) && (ETA_R() <= 1.0)))
    {
        throw std::invalid_argument("eta_r must be between 0.0 and 1.0. You entered: " + std::to_string(ETA_R()));
    }
    if (Spacing() == std::string{})
    {
        throw std::invalid_argument("BSpline spacing must be non-empty. You entered: " + spacing);
    }
}

void BSpline::buildKnots(const Box& box)
{
    if (spacing == "linear")
    {
        buildLinearKnots(box);
    }
    else
    {
        throw std::invalid_argument("No matching BSpline Spacing. You entered: " + spacing);
    }


    buildR0();
    buildComplexKnots();
}

void BSpline::buildComplexKnots()
{
    int N_knots = n_bspline + order;

    complex_knots.resize(N_knots);

    for (int idx = 0; idx < N_knots; ++idx) 
    {
        complex_knots[idx] = ecs_x(knots[idx]);
    }
}

void BSpline::buildLinearKnots(const Box& box)
{
    int N_knots   = n_bspline + order;
    int N_middle  = N_knots - 2 * order;
    double step   = box.GridSize() / (N_middle - 1);

    knots.resize(N_knots);
    int start_mid = order;
    int end_mid   = N_knots - order - 1;

    for (int idx = 0; idx < N_knots; ++idx) 
    {
        if (idx < start_mid) 
        {
            knots[idx] = 0.0;
        }
        else if (idx > end_mid) 
        {
            knots[idx] = box.GridSize();
        }
        else 
        {
            int j = idx - start_mid;
            knots[idx] = j * step;
        }
    }
}

std::complex<double> BSpline::ecs_x(double x) const
{
    if (x < R0) 
    {
        return std::complex<double>{x, 0.0};
    }
    else 
    {
        return R0 + (x - R0) * std::exp(std::complex<double>{0, eta});
    }
}

std::complex<double> BSpline::ecs_w(double x, double w) const
{
    if (x < R0) 
    {
        return std::complex<double>{w, 0.0};
    }
    else 
    {
        return w * std::exp(std::complex<double>{0, eta});
    }
}

void BSpline::buildR0()
{
    double target   = R0_r * knots.back();
    double min_val  = std::abs(knots[0] - target);
    double knot_val = knots[0];

    for (size_t idx = 1; idx < knots.size(); ++idx) {
        double diff = std::abs(knots[idx] - target);
        if (diff < min_val) {
            min_val  = diff;
            knot_val = knots[idx];
        }
    }

    R0 = knot_val;
}

std::complex<double> BSpline::B(int degree, int i, std::complex<double> x) const
{
    if (degree == 0)
    {
        return (complex_knots[i].real() <= x.real() && x.real() < complex_knots[i + 1].real() ? 1.0 : 0.0);
    }

    std::complex<double> denom1 = complex_knots[i + degree] - complex_knots[i];
    std::complex<double> denom2 = complex_knots[i + degree + 1] - complex_knots[i + 1];

    std::complex<double> term1 = 0.0;
    std::complex<double> term2 = 0.0;

    if (denom1.real() > 0)
    {
        term1 = (x - complex_knots[i]) / denom1 * B(degree - 1, i, x);
    }
    if (denom2.real() > 0)
    {
        term2 = (complex_knots[i + degree + 1] - x) / denom2 * B(degree - 1, i + 1, x);
    }

    return term1 + term2;
}

std::complex<double> BSpline::dB( int degree, int i, std::complex<double> x) const
{
    if (degree == 0)
    {
        return 0.0;
    }

    std::complex<double> denom1 = complex_knots[i + degree] - complex_knots[i];
    std::complex<double> denom2 = complex_knots[i + degree + 1] - complex_knots[i + 1];

    std::complex<double> term1 = 0.0;
    std::complex<double> term2 = 0.0;

    if (denom1.real() > 0)
    {
        term1 = std::complex<double>(degree) / denom1 * B(degree - 1, i, x);
    }
    if (denom2.real() > 0)
    {
        term2 = -std::complex<double>(degree) / denom2 * B(degree - 1, i + 1, x);
    }

    return term1 + term2;
}

std::complex<double> BSpline::integrateMatrixElement(int i, int j,std::function<std::complex<double>(int, int, std::complex<double>)> integrand,bool use_ecs) const
{
    std::complex<double> total{0.0,0.0};

    int lower = std::min(i, j);
    int upper = std::max(i, j);

    for (int k = lower; k <= upper + degree; ++k)
    {
        double a = knots[k];
        double b = knots[k + 1];


        if (a == b)
        {
            continue;
        }

        double half_b_minus_a = 0.5 * (b - a);
        double half_b_plus_a = 0.5 * (b + a);


        for (size_t r = 0; r < roots.size(); ++r)
        {
            double x_val = half_b_minus_a * roots[r] + half_b_plus_a;
            double weight_val = weights[r];

            if (use_ecs)
            {
                std::complex<double> x = ecs_x(x_val);
                std::complex<double> weight = ecs_w(x_val, weight_val) * half_b_minus_a;
                std::complex<double> integrand_val = integrand(i, j, x);
                total += weight * integrand_val;
            }
            else
            {
                std::complex<double> x = x_val;
                std::complex<double> weight = weight_val* half_b_minus_a;
                std::complex<double> integrand_val = integrand(i, j, x);
                total += weight * integrand_val;
            }
        }
    }

    return total;
}

void BSpline::dumpTo(const Box& box, const std::string& directory, int rank)
{
    if (rank == 0) 
    {
        std::string filename = directory + "/bsplines.txt";

        std::ofstream outFile(filename);

        if (!outFile) 
        {
            std::cerr << "Error opening file: " << filename << '\n';
        }

        for (int spline = 0; spline < n_bspline; ++spline) 
        {
            for (int ridx = 0; ridx < box.Nr(); ++ridx) 
            {   
                std::complex<double> value = B(spline, ecs_x(ridx * box.GridSpacing()));
                outFile << value.real() << " " << value.imag() << '\n';
            }
        }
        outFile.close();

        std::string filename2 = directory + "/dbsplines.txt";

        std::ofstream outFile2(filename2);

        if (!outFile2) 
        {
            std::cerr << "Error opening file: " << filename2 << '\n';
        }

        for (int spline = 0; spline < n_bspline; ++spline) 
        {
            for (int ridx = 0; ridx < box.Nr(); ++ridx) 
            {
                std::complex<double> value = dB(spline, ecs_x(ridx * box.GridSpacing()));
                outFile2 << value.real() << " " << value.imag() << '\n';
            }
        }
        outFile2.close();

        std::string filename3 = directory + "/bspline_metadata.txt";

        std::ofstream outFile3(filename3);

        if (!outFile3) 
        {
            std::cerr << "Error opening file: " << filename2 << '\n';
        }

        outFile3 << box.Nr() << " " << box.GridSpacing() << '\n';

        outFile3.close();
    }
}
