
#include "common.h"

#include "Basis.h"
#include "Box.h"
#include "PetscWrappers/PetscMat.h"
#include "BSplines.h"

const std::unordered_map<int, std::pair<std::vector<double>, std::vector<double>>> Basis::GaussLegendre::RootsAndWeights = 
{
    {2, {{-0.57735027, 0.57735027}, {1, 1}}},
    {3, {{-0.77459667, 0.0, 0.77459667}, {0.55555556, 0.88888889, 0.55555556}}},
    {4, {{-0.86113631, -0.33998104, 0.33998104, 0.86113631}, {0.34785485, 0.65214515, 0.65214515, 0.34785485}}},
    {5, {{-0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985}, {0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689}}},
    {6, {{-0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951}, {0.17132449, 0.36076157, 0.46791393, 0.46791393, 0.36076157, 0.17132449}}},
    {7, {{-0.94910791, -0.74153119, -0.40584515, 0, 0.40584515, 0.74153119, 0.94910791}, {0.12948497, 0.27970539, 0.38183005, 0.41795918, 0.38183005, 0.27970539, 0.12948497}}},
    {8, {{-0.96028986, -0.79666648, -0.52553241, -0.18343464, 0.18343464, 0.52553241, 0.79666648, 0.96028986}, {0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854}}}
};

void Basis::buildParameters()
{
    m_degree = getOrder() - 1;
    m_roots = GaussLegendre::RootsAndWeights.at(getOrder()).first;
    m_weights = GaussLegendre::RootsAndWeights.at(getOrder()).second;
            
    m_eta *= M_PI;
}

void Basis::buildKnots(const Box& box)
{
    if (getSpacing() == "linear")
    {
        buildLinearKnots(box);
    }
   
    buildR0();
    buildComplexKnots();
}

void Basis::buildComplexKnots()
{
    int N_knots = getNbasis() + getOrder();

    m_complexKnots.resize(N_knots);

    for (int idx = 0; idx < N_knots; ++idx) 
    {
        m_complexKnots[idx] = ecs_x(std::real(m_knots[idx]));
    }
}

void Basis::buildLinearKnots(const Box& box)
{
    int N_knots   = getNbasis() + getOrder();
    int leftMult  = getOrder() - 2;
    int rightMult = getOrder() - 2;
    int N_middle  = N_knots - leftMult - rightMult;   
    double step   = box.getGridSize() / (N_middle - 1);

    m_knots.clear();
    m_knots.reserve(N_knots);

    for (int i = 0; i < leftMult; ++i)
    {
        m_knots.push_back(0.0);
    }
        
    for (int j = 0; j < N_middle; ++j)
    {
        m_knots.push_back(j * step);
    }
        
    for (int i = 0; i < rightMult; ++i)
    {
        m_knots.push_back(box.getGridSize());
    }

}

void Basis::buildR0()
{
    double target   = getR0() * std::real(getKnots().back());
    double min_val  = std::abs(getKnots()[0] - target);
    double knot_val = std::real(getKnots()[0]);

    for (size_t idx = 1; idx < getKnots().size(); ++idx) {
        double diff = std::abs(getKnots()[idx] - target);
        if (diff < min_val) {
            min_val  = diff;
            knot_val = std::real(getKnots()[idx]);
        }
    }

    m_R0 = knot_val;
}

std::complex<double> Basis::integrateMatrixElement(int i, int j,MatrixIntegrand integrand,bool use_ecs) const
{
    std::complex<double> total{0.0,0.0};

    int lower = std::min(i, j);
    int upper = std::max(i, j);

    for (int k = lower; k <= upper + getDegree(); ++k)
    {
        double a = std::real(getKnots()[k]);
        double b = std::real(getKnots()[k + 1]);


        if (a == b)
        {
            continue;
        }

        double half_b_minus_a = 0.5 * (b - a);
        double half_b_plus_a = 0.5 * (b + a);


        for (size_t r = 0; r < getRoots().size(); ++r)
        {
            double x_val = half_b_minus_a * getRoots()[r] + half_b_plus_a;
            double weight_val = getWeights()[r];

            if (use_ecs)
            {
                std::complex<double> x = ecs_x(x_val);
                std::complex<double> weight = ecs_w(x_val, weight_val) * half_b_minus_a;
                std::complex<double> integrand_val = (*integrand)(i, j,getDegree(), x,getComplexKnots());
                total += weight * integrand_val;
            }
            else
            {
                std::complex<double> x = x_val;
                std::complex<double> weight = weight_val* half_b_minus_a;
                std::complex<double> integrand_val = (*integrand)(i, j,getDegree(), x,getKnots());
                total += weight * integrand_val;
            }
        }
    }

    return total;
}


// void Basis::dumpTo(const Box& box, const std::string& directory, int rank)
// {
//     if (rank == 0) 
//     {
//         std::string filename = directory + "/Basiss.txt";

//         std::ofstream outFile(filename);

//         if (!outFile) 
//         {
//             std::cerr << "Error opening file: " << filename << '\n';
//         }

//         for (int spline = 0; spline < getNbasis(); ++spline) 
//         {
//             for (int ridx = 0; ridx < box.getNr(); ++ridx) 
//             {   
//                 std::complex<double> value = BSplines::B(spline,getDegree(), ecs_x(ridx * box.getGridSpacing()),complex_knots);
//                 outFile << value.real() << " " << value.imag() << '\n';
//             }
//         }
//         outFile.close();

//         std::string filename2 = directory + "/dBasiss.txt";

//         std::ofstream outFile2(filename2);

//         if (!outFile2) 
//         {
//             std::cerr << "Error opening file: " << filename2 << '\n';
//         }

//         for (int spline = 0; spline < getNbasis(); ++spline) 
//         {
//             for (int ridx = 0; ridx < box.getNr(); ++ridx) 
//             {
//                 std::complex<double> value = BSplines::dB(spline,getDegree(), ecs_x(ridx * box.getGridSpacing()),complex_knots);
//                 outFile2 << value.real() << " " << value.imag() << '\n';
//             }
//         }
//         outFile2.close();

//         std::string filename3 = directory + "/Basis_metadata.txt";

//         std::ofstream outFile3(filename3);

//         if (!outFile3) 
//         {
//             std::cerr << "Error opening file: " << filename2 << '\n';
//         }

//         outFile3 << box.getNr() << " " << box.getGridSpacing() << '\n';

//         outFile3.close();
//     }
// }












