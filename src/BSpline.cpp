
#include "BSpline.h"
#include <complex>
#include <fstream>
#include <iostream>
#include "Box.h"

void BSpline::buildLinearKnots(const Box& box)
{
    int N_knots   = n_bspline + order;
    int N_middle  = N_knots - 2 * order;
    // use (N_middle - 1) here to get exactly N_middle interior knots spanning [0,GridSize]
    double step   = box.GridSize() / (N_middle - 1);

    knots.resize(N_knots);
    int start_mid = order;
    int end_mid   = N_knots - order - 1;

    for (int idx = 0; idx < N_knots; ++idx) {
        if (idx < start_mid) {
            knots[idx] = 0.0;
        }
        else if (idx > end_mid) {
            knots[idx] = box.GridSize();
        }
        else {
            int j = idx - start_mid;
            knots[idx] = j * step;
        }
    }

    buildR0();
    buildComplexKnots();
}

void BSpline::buildComplexKnots()
{
    int N_knots = n_bspline + order;

    complex_knots.clear();
    complex_knots.reserve(N_knots);

    for (int idx = 0; idx < N_knots; ++idx) {
        complex_knots.push_back(ecs_x(knots[idx]));
    }
}

std::complex<double> BSpline::ecs_x(double x) const
{
    if (x < R0) {
        return std::complex<double>(x, 0.0);
    }
    else {
        return R0 + (x - R0) * std::exp(std::complex<double>(0, M_PI * eta));
    }
}

std::complex<double> BSpline::ecs_w(double x, double w) const
{
    if (x < R0) {
        return std::complex<double>(w, 0.0);
    }
    else {
        return w * std::exp(std::complex<double>(0, M_PI * eta));
    }
}

void BSpline::buildR0()
{
    double target   = R0 * knots.back();
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

std::complex<double> BSpline::B(int i, std::complex<double> x) const
{
    std::vector<std::complex<double>> N(degree + 1);

    // zeroth-order using complex_knots.real() for support test
    for (int j = 0; j <= degree; ++j) {
        double t0 = complex_knots[i + j].real();
        double t1 = complex_knots[i + j + 1].real();
        N[j] = (t0 <= x.real() && x.real() < t1)
               ? std::complex<double>(1.0, 0.0)
               : std::complex<double>(0.0, 0.0);
    }

    for (int d = 1; d <= degree; ++d) {
        for (int j = 0; j <= degree - d; ++j) {
            auto left_den  = complex_knots[i + j + d]     - complex_knots[i + j];
            auto right_den = complex_knots[i + j + d + 1] - complex_knots[i + j + 1];

            std::complex<double> t1(0.0, 0.0), t2(0.0, 0.0);

            if (std::abs(left_den)  > 0.0)
                t1 = (x - complex_knots[i + j])       / left_den  * N[j];
            if (std::abs(right_den) > 0.0)
                t2 = (complex_knots[i + j + d + 1] - x) / right_den * N[j + 1];

            N[j] = t1 + t2;
        }
    }

    return N[0];
}

std::complex<double> BSpline::dB(int i, std::complex<double> x) const
{
    const int p = degree;
    if (p == 0)
        return std::complex<double>(0.0, 0.0);

    auto compute_basis = [&](int ii) {
        std::vector<std::complex<double>> N(p);
        // zeroth-order using complex_knots.real() for support test
        for (int j = 0; j < p; ++j) {
            double t0 = complex_knots[ii + j].real();
            double t1 = complex_knots[ii + j + 1].real();
            N[j] = (t0 <= x.real() && x.real() < t1)
                   ? std::complex<double>(1.0, 0.0)
                   : std::complex<double>(0.0, 0.0);
        }
        // build up to p-1
        for (int d = 1; d < p; ++d) {
            for (int j = 0; j < p - d; ++j) {
                auto ld = complex_knots[ii + j + d]     - complex_knots[ii + j];
                auto rd = complex_knots[ii + j + d + 1] - complex_knots[ii + j + 1];
                std::complex<double> a(0.0, 0.0), b(0.0, 0.0);
                if (std::abs(ld) > 0.0)
                    a = (x - complex_knots[ii + j])       / ld * N[j];
                if (std::abs(rd) > 0.0)
                    b = (complex_knots[ii + j + d + 1] - x) / rd * N[j + 1];
                N[j] = a + b;
            }
        }
        return N[0];
    };

    std::complex<double> N_i   = compute_basis(i);
    std::complex<double> N_ip1 = compute_basis(i + 1);

    auto denom1 = complex_knots[i + p]     - complex_knots[i];
    auto denom2 = complex_knots[i + p + 1] - complex_knots[i + 1];
    double pd = static_cast<double>(p);

    std::complex<double> term1(0.0, 0.0), term2(0.0, 0.0);
    if (std::abs(denom1) > 0.0)
        term1 = pd / denom1 * N_i;
    if (std::abs(denom2) > 0.0)
        term2 = pd / denom2 * N_ip1;

    return term1 - term2;
}

void BSpline::dumpTo(const Box& box, const std::string& directory, int rank)
{
    if (rank == 0) {
        std::string filename = directory + "/bsplines.txt";
        std::ofstream outFile(filename);
        if (!outFile) {
            std::cerr << "Error opening file: " << filename << '\n';
        }
        for (int spline = 0; spline < n_bspline; ++spline) {
            for (int ridx = 0; ridx < box.Nr(); ++ridx) {
                outFile << B(spline, ecs_x(ridx * box.GridSpacing()));
            }
            outFile << '\n';
        }
        outFile.close();

        std::string filename2 = directory + "/dbsplines.txt";
        std::ofstream outFile2(filename2);
        if (!outFile2) {
            std::cerr << "Error opening file: " << filename2 << '\n';
        }
        for (int spline = 0; spline < n_bspline; ++spline) {
            for (int ridx = 0; ridx < box.Nr(); ++ridx) {
                outFile2 << dB(spline, ecs_x(ridx * box.GridSpacing()));
            }
            outFile2 << '\n';
        }
        outFile2.close();
    }
}
