#pragma once
#include "nlohmann/json.hpp"
#include "Box.h"
#include <complex>
#include "Potentials.h"
#include "GaussLegendre.h"
#include "PetscWrappers/PetscMat.h"

class BSpline
{
    public:
        BSpline() = delete;
        explicit BSpline(const nlohmann::json& input_file, const Box& box)
        : n_bspline{input_file.at("BSpline").at("n_bspline")}
        , order   {input_file.at("BSpline").at("order")}
        , R0_r      {input_file.at("BSpline").at("R0_r")}
        , eta_r     {input_file.at("BSpline").at("eta_r")}
        , spacing {input_file.at("BSpline").at("spacing")}
        {
            validateInput();

            degree = order - 1;
            roots = GS::Gauss.at(order).first;
            weights = GS::Gauss.at(order).second;
            
            eta = eta_r * M_PI;

            
            buildKnots(box);
        }

        int NBSpline() const {return n_bspline;}
        int Order() const {return order;}
        int Degree() const {return degree;}
        double R0_R() const {return R0_r;}
        double ETA_R() const {return eta_r;}
        const std::string& Spacing() const {return spacing;}
        const std::vector<double> Knots() const {return knots;}
        const std::vector<std::complex<double>> ComplexKnots() const {return complex_knots;}

        std::complex<double> BTest(int i, std::complex<double> x) const;
        std::complex<double> dBTest(int i, std::complex<double> x) const;
        std::complex<double> B(int degree, int i, std::complex<double> x) const;
        std::complex<double> B(int i, std::complex<double> x) const {return B(degree, i, x);}
        std::complex<double> dB(int degree, int i, std::complex<double> x) const;
        std::complex<double> dB(int i, std::complex<double> x) const {return dB(degree, i ,x);}
        void dumpTo(const Box& box, const std::string& directory, int rank);

    private:
        // Member List Initialized
        int n_bspline{};
        int order{};
        int degree{};
        double R0_r{};
        double eta_r{};
        std::string spacing{};

        // Default Initialized
        double R0;
        double eta;
        std::vector<double> knots;
        std::vector<std::complex<double>> complex_knots;
        std::vector<double> roots;
        std::vector<double> weights;

        // Member Functions
        void validateInput();
    public:
        std::complex<double> ecs_x(double x) const;
        std::complex<double> ecs_w(double x, double w) const;
    private:
        void buildKnots(const Box& box);
        void buildLinearKnots(const Box& box);
        void buildComplexKnots();
        void buildR0();

        std::complex<double> overlapIntegrand(int i, int j, std::complex<double> x) const {return B(i, x) * B(j, x);}
        std::complex<double> kineticIntegrand(int i, int j, std::complex<double> x) const {return 0.5 * dB(i, x) * dB(j, x);}
        std::complex<double> invrIntegrand(int i, int j, std::complex<double> x) const {return B(i, x) * B(j, x) / (x + 1E-25);}
        std::complex<double> invr2Integrand(int i, int j, std::complex<double> x) const {return B(i, x) * B(j, x) / (x*x + 1E-25);}
        std::complex<double> derIntegrand(int i, int j, std::complex<double> x) const {return B(i, x) * dB(j,x);}
        std::complex<double> HIntegrand(int i, int j, std::complex<double> x) const {return B(i, x) * B(j, x) * Potentials::hydrogenPotential(x);}
        std::complex<double> integrateMatrixElement(int i, int j,std::function<std::complex<double>(int, int, std::complex<double>)> integrand,bool use_ecs) const;

        void PopulateMatrix(Matrix& matrix,std::function<std::complex<double>(int, int, std::complex<double>)> integrand,bool use_ecs);
};






    
