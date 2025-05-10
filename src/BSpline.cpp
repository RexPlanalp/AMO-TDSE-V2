
#include "common.h"

#include "BSpline.h"
#include "Box.h"


const std::unordered_map<int, std::pair<std::vector<double>, std::vector<double>>> BSpline::GaussLegendre::RootsAndWeights = 
{
    {2, {{-0.57735027, 0.57735027}, {1, 1}}},
    {3, {{-0.77459667, 0.0, 0.77459667}, {0.55555556, 0.88888889, 0.55555556}}},
    {4, {{-0.86113631, -0.33998104, 0.33998104, 0.86113631}, {0.34785485, 0.65214515, 0.65214515, 0.34785485}}},
    {5, {{-0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985}, {0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689}}},
    {6, {{-0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951}, {0.17132449, 0.36076157, 0.46791393, 0.46791393, 0.36076157, 0.17132449}}},
    {7, {{-0.94910791, -0.74153119, -0.40584515, 0, 0.40584515, 0.74153119, 0.94910791}, {0.12948497, 0.27970539, 0.38183005, 0.41795918, 0.38183005, 0.27970539, 0.12948497}}},
    {8, {{-0.96028986, -0.79666648, -0.52553241, -0.18343464, 0.18343464, 0.52553241, 0.79666648, 0.96028986}, {0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854}}}
};

void BSpline::buildKnots(const Box& box)
{
    if (spacing == "linear")
    {
        buildLinearKnots(box);
    }
   
    buildR0();
    buildComplexKnots();
}

void BSpline::buildComplexKnots()
{
    int N_knots = getNbasis() + getOrder();

    complex_knots.resize(N_knots);

    for (int idx = 0; idx < N_knots; ++idx) 
    {
        complex_knots[idx] = ecs_x(knots[idx]);
    }
}

void BSpline::buildLinearKnots(const Box& box)
{
    int N_knots   = getNbasis() + getOrder();
    int N_middle  = N_knots - 2 * getOrder();
    double step   = box.getGridSize() / (N_middle - 1);

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
            knots[idx] = box.getGridSize();
        }
        else 
        {
            int j = idx - start_mid;
            knots[idx] = j * step;
        }
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

std::complex<double> BSpline::integrateMatrixElement(int i, int j, BSpline::MatrixIntegrand integrand,bool use_ecs) const
{
    std::complex<double> total{0.0,0.0};

    int lower = std::min(i, j);
    int upper = std::max(i, j);

    for (int k = lower; k <= upper + getDegree(); ++k)
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
                std::complex<double> integrand_val = (this->*integrand)(i, j, x);
                total += weight * integrand_val;
            }
            else
            {
                std::complex<double> x = x_val;
                std::complex<double> weight = weight_val* half_b_minus_a;
                std::complex<double> integrand_val = (this->*integrand)(i, j, x);
                total += weight * integrand_val;
            }
        }
    }

    return total;
}

void BSpline::printConfiguration(int rank)
{
    if (rank == 0)
    {
        std::cout << "BSpline Configuration: " << "\n\n";
        std::cout << "nbasis: " << getNbasis() <<  "\n\n";
        std::cout << "order: " << getOrder() <<  "\n\n";
        std::cout << "degree: " << getDegree() <<  "\n\n";
        std::cout << "spacing: " << getSpacing() <<  "\n\n";
    }
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

        for (int spline = 0; spline < getNbasis(); ++spline) 
        {
            for (int ridx = 0; ridx < box.getNr(); ++ridx) 
            {   
                std::complex<double> value = B(spline, ecs_x(ridx * box.getGridSpacing()));
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

        for (int spline = 0; spline < getNbasis(); ++spline) 
        {
            for (int ridx = 0; ridx < box.getNr(); ++ridx) 
            {
                std::complex<double> value = dB(spline, ecs_x(ridx * box.getGridSpacing()));
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

        outFile3 << box.getNr() << " " << box.getGridSpacing() << '\n';

        outFile3.close();
    }
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

// Matrix BSpline::PopulateMatrix(std::function<std::complex<double>(int, int, std::complex<double>)> integrand,bool use_ecs) const
// {   
//     Matrix matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,NBSpline(),NBSpline(),2*Degree() + 1};

   

//     for (int i = matrix.getStart(); i < matrix.getEnd(); i++) 
//     {
//         int col_start = std::max(0, i - order + 1);
//         int col_end = std::min(NBSpline(), i + order);

//         for (int j = col_start; j < col_end; j++) 
//         {
//             std::complex<double> result = integrateMatrixElement(i, j,integrand,use_ecs);
//             MatSetValue(matrix.get(), i, j, result, INSERT_VALUES); 
//         }
//     }
//     matrix.assemble();
//     return matrix;
// }













// std::complex<double> BSpline::B(int degree, int i, std::complex<double> x) const
// {
//     if (degree == 0)
//     {
//         return (complex_knots[i].real() <= x.real() && x.real() < complex_knots[i + 1].real() ? 1.0 : 0.0);
//     }

//     std::complex<double> denom1 = complex_knots[i + degree] - complex_knots[i];
//     std::complex<double> denom2 = complex_knots[i + degree + 1] - complex_knots[i + 1];

//     std::complex<double> term1 = 0.0;
//     std::complex<double> term2 = 0.0;

//     if (denom1.real() > 0)
//     {
//         term1 = (x - complex_knots[i]) / denom1 * B(degree - 1, i, x);
//     }
//     if (denom2.real() > 0)
//     {
//         term2 = (complex_knots[i + degree + 1] - x) / denom2 * B(degree - 1, i + 1, x);
//     }

//     return term1 + term2;
// }

// std::complex<double> BSpline::dB( int degree, int i, std::complex<double> x) const
// {
//     if (degree == 0)
//     {
//         return 0.0;
//     }

//     std::complex<double> denom1 = complex_knots[i + degree] - complex_knots[i];
//     std::complex<double> denom2 = complex_knots[i + degree + 1] - complex_knots[i + 1];

//     std::complex<double> term1 = 0.0;
//     std::complex<double> term2 = 0.0;

//     if (denom1.real() > 0)
//     {
//         term1 = std::complex<double>(degree) / denom1 * B(degree - 1, i, x);
//     }
//     if (denom2.real() > 0)
//     {
//         term2 = -std::complex<double>(degree) / denom2 * B(degree - 1, i + 1, x);
//     }

//     return term1 + term2;
// }