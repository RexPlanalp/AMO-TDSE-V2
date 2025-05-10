#pragma once

#include "Input.h"

#include "Box.h"


class BSpline
{   
    using MatrixIntegrand = std::complex<double> (BSpline::*)(int, int, std::complex<double>) const;
    struct GaussLegendre
    {
        static const std::unordered_map<int, std::pair<std::vector<double>, std::vector<double>>> RootsAndWeights;
    };



    public:
        explicit BSpline(const Input& input)
        : nbasis{input.getJSON().at("BSpline").at("nbasis")}
        , order{input.getJSON().at("BSpline").at("order")}
        , R0{input.getJSON().at("BSpline").at("R0r")}
        , eta{input.getJSON().at("BSpline").at("etar")}
        , spacing{input.getJSON().at("BSpline").at("spacing")}
        {
            degree = order - 1;
            roots = GaussLegendre::RootsAndWeights.at(order).first;
            weights = GaussLegendre::RootsAndWeights.at(order).second;
            
            eta *= M_PI;
        }

        int getNbasis() const {return nbasis;}
        int getOrder() const {return order;}
        int getDegree() const {return degree;}
        double getR0() const {return R0;}
        double getEta() const {return eta;}
        const std::string& getSpacing() const {return spacing;}
        const std::vector<double> getKnots() const {return knots;}
        const std::vector<std::complex<double>> getComplexKnots() const {return complex_knots;}

        std::complex<double> B(int i, std::complex<double> x) const;
        std::complex<double> dB(int i, std::complex<double> x) const;

        std::complex<double> integrateMatrixElement(int i, int j, MatrixIntegrand integrand,bool use_ecs) const;

        void printConfiguration(int rank);
        void dumpTo(const Box& box, const std::string& directory, int rank);
        void buildKnots(const Box& box);
        

        // std::complex<double> B(int degree, int i, std::complex<double> x) const;
        // std::complex<double> B(int i, std::complex<double> x) const {return B(degree, i, x);}
        // std::complex<double> dB(int degree, int i, std::complex<double> x) const;
        // std::complex<double> dB(int i, std::complex<double> x) const {return dB(degree, i ,x);}

    private:
        // Member List Initialized
        int nbasis{};
        int order{};
        double R0{};
        double eta{};
        std::string spacing{};

        // Derived
        int degree{};
        std::vector<double> knots{};
        std::vector<std::complex<double>> complex_knots{};
        std::vector<double> roots{};
        std::vector<double> weights{};

        // Member Functions
        std::complex<double> ecs_x(double x) const {return (x < R0) ? (std::complex<double>{x, 0.0}) : R0 + (x - R0) * std::exp(std::complex<double>{0, eta});}
        std::complex<double> ecs_w(double x, double w) const { return (x < R0) ? (std::complex<double>{w, 0.0}) : w * std::exp(std::complex<double>{0, eta});}
        void buildLinearKnots(const Box& box);
        void buildComplexKnots();
        void buildR0();

        std::complex<double> overlapIntegrand(int i, int j, std::complex<double> x) const {return B(i, x) * B(j, x);}
        std::complex<double> kineticIntegrand(int i, int j, std::complex<double> x) const {return 0.5 * dB(i, x) * dB(j, x);}
        std::complex<double> invrIntegrand(int i, int j, std::complex<double> x) const {return B(i, x) * B(j, x) / (x + 1E-25);}
        std::complex<double> invr2Integrand(int i, int j, std::complex<double> x) const {return B(i, x) * B(j, x) / (x*x + 1E-25);}
        std::complex<double> derIntegrand(int i, int j, std::complex<double> x) const {return B(i, x) * dB(j,x);}



        // std::complex<double> HIntegrand(int i, int j, std::complex<double> x) const {return B(i, x) * B(j, x) * Potentials::hydrogenPotential(x);}
        

        // Matrix PopulateMatrix(std::function<std::complex<double>(int, int, std::complex<double>)> integrand,bool use_ecs) const;


        
};






    
