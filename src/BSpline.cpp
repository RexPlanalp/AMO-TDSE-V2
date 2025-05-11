
#include "common.h"

#include "BSpline.h"
#include "Box.h"
#include "PetscWrappers/PetscMat.h"

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
        complex_knots[idx] = ecs_x(std::real(knots[idx]));
    }
}



void BSpline::buildLinearKnots(const Box& box)
{
    int N_knots   = getNbasis() + getOrder();
    int leftMult  = getOrder() - 2;
    int rightMult = getOrder() - 2;
    int N_middle  = N_knots - leftMult - rightMult;   
    double step   = box.getGridSize() / (N_middle - 1);

    knots.clear();
    knots.reserve(N_knots);

    for (int i = 0; i < leftMult; ++i)
    {
        knots.push_back(0.0);
    }
        
    for (int j = 0; j < N_middle; ++j)
    {
        knots.push_back(j * step);
    }
        
    for (int i = 0; i < rightMult; ++i)
    {
        knots.push_back(box.getGridSize());
    }

}

void BSpline::buildR0()
{
    double target   = R0 * std::real(knots.back());
    double min_val  = std::abs(knots[0] - target);
    double knot_val = std::real(knots[0]);

    for (size_t idx = 1; idx < knots.size(); ++idx) {
        double diff = std::abs(knots[idx] - target);
        if (diff < min_val) {
            min_val  = diff;
            knot_val = std::real(knots[idx]);
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
        double a = std::real(knots[k]);
        double b = std::real(knots[k + 1]);


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
                std::complex<double> integrand_val = (this->*integrand)(i, j, x,complex_knots);
                total += weight * integrand_val;
            }
            else
            {
                std::complex<double> x = x_val;
                std::complex<double> weight = weight_val* half_b_minus_a;
                std::complex<double> integrand_val = (this->*integrand)(i, j, x,knots);
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
                std::complex<double> value = B(spline, ecs_x(ridx * box.getGridSpacing()),complex_knots);
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
                std::complex<double> value = dB(spline, ecs_x(ridx * box.getGridSpacing()),complex_knots);
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

std::complex<double> BSpline::dB(int i, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) const
{
    const int p = degree;
    if (p == 0)
        return std::complex<double>(0.0, 0.0);

    auto compute_basis = [&](int ii) {
        std::vector<std::complex<double>> N(p);
        // zeroth-order using complex_knots.real() for support test
        for (int j = 0; j < p; ++j) {
            double t0 = std::real(localKnots[ii + j]);
            double t1 = std::real(localKnots[ii + j + 1]);
            N[j] = (t0 <= x.real() && x.real() < t1)
                   ? std::complex<double>(1.0, 0.0)
                   : std::complex<double>(0.0, 0.0);
        }
        // build up to p-1
        for (int d = 1; d < p; ++d) {
            for (int j = 0; j < p - d; ++j) {
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

std::complex<double> BSpline::B(int i, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) const
{
    std::vector<std::complex<double>> N(degree + 1);

    // zeroth-order using complex_knots.real() for support test
    for (int j = 0; j <= degree; ++j) {
        double t0 = std::real(localKnots[i + j]);
        double t1 = std::real(localKnots[i + j + 1]);
        N[j] = (t0 <= std::real(x) && std::real(x) < t1)
               ? std::complex<double>(1.0, 0.0)
               : std::complex<double>(0.0, 0.0);
    }

    for (int d = 1; d <= degree; ++d) {
        for (int j = 0; j <= degree - d; ++j) {
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

Matrix BSpline::PopulateMatrix(BSpline::MatrixIntegrand integrand,bool use_ecs) const
{   
    Matrix matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getNbasis(),getNbasis(),2*getDegree() + 1};

   

    for (int i = matrix.getStart(); i < matrix.getEnd(); i++) 
    {
        int col_start = std::max(0, i - getOrder() + 1);
        int col_end = std::min(getNbasis(), i + getOrder());

        for (int j = col_start; j < col_end; j++) 
        {
            std::complex<double> result = integrateMatrixElement(i, j,integrand,use_ecs);
            MatSetValue(matrix.get(), i, j, result, INSERT_VALUES); 
        }
    }
    matrix.assemble();
    return matrix;
}











