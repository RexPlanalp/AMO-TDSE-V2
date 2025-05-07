#pragma once
#include "nlohmann/json.hpp"
#include "Box.h"
#include <complex>
#include "Potentials.h"
#include "GaussLegendre.h"

class BSpline
{
    
        
       


    public:
        BSpline() = delete;
        explicit BSpline(const nlohmann::json& input_file, const Box& box)
        : n_bspline{input_file.at("BSpline").at("n_bspline")}
        , order   {input_file.at("BSpline").at("order")}
        , R0      {input_file.at("BSpline").at("R0_r")}
        , eta     {input_file.at("BSpline").at("eta_r")}
        , spacing {input_file.at("BSpline").at("spacing")}
        {
            degree = order - 1;
            roots = GS::Gauss.at(order).first;
            weights = GS::Gauss.at(order).second;
            
            buildKnots(box);
        }

        void buildKnots(const Box& box);
        

        std::complex<double> ecs_x(double x) const;
        std::complex<double> ecs_w(double x, double w) const;
        std::complex<double> B(int i, std::complex<double> x) const;
        std::complex<double> dB(int i, std::complex<double> x) const;

        void dumpTo(const Box& box, const std::string& directory, int rank);

    private:
        
        // Member List Initialized
        int n_bspline{};
        int order{};
        int degree{};
        double R0{};
        double eta{};
        std::string spacing{};

        // Default Initialized
        std::vector<double> knots;
        std::vector<std::complex<double>> complex_knots;
        std::vector<double> roots;
        std::vector<double> weights;

        // Member Functions
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
};




    
