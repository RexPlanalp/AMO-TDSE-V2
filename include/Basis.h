#pragma once

#include "Input.h"

#include "Box.h"
#include "PetscWrappers/PetscMat.h"

#include "Potentials.h"

using MatrixIntegrand = std::complex<double> (*)(int, int,int, std::complex<double>, const std::vector<std::complex<double>>&);

class Basis
{   
    
    
    struct GaussLegendre
    {
        static const std::unordered_map<int, std::pair<std::vector<double>, std::vector<double>>> RootsAndWeights;
    };

    public:
        explicit Basis(const Input& input)
        : nbasis{input.getJSON().at("Basis").at("nbasis")}
        , order{input.getJSON().at("Basis").at("order")}
        , R0{input.getJSON().at("Basis").at("R0r")}
        , eta{input.getJSON().at("Basis").at("etar")}
        , spacing{input.getJSON().at("Basis").at("spacing")}
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
        const std::vector<std::complex<double>> getKnots() const {return knots;}
        const std::vector<std::complex<double>> getComplexKnots() const {return complex_knots;}
        
        std::complex<double> Basis::integrateMatrixElement(int i, int j,MatrixIntegrand integrand,bool use_ecs) const;

        void printConfiguration(int rank);
        void dumpTo(const Box& box, const std::string& directory, int rank);
        void buildKnots(const Box& box);


        
    

    private:
        // Member List Initialized
        int nbasis{};
        int order{};
        double R0{};
        double eta{};
        std::string spacing{};

        // Derived
        int degree{};
        std::vector<std::complex<double>> knots{};
        std::vector<std::complex<double>> complex_knots{};
        std::vector<double> roots{};
        std::vector<double> weights{};

        // Member Functions
        std::complex<double> ecs_x(double x) const {return (x < R0) ? (std::complex<double>{x, 0.0}) : R0 + (x - R0) * std::exp(std::complex<double>{0, eta});}
        std::complex<double> ecs_w(double x, double w) const { return (x < R0) ? (std::complex<double>{w, 0.0}) : w * std::exp(std::complex<double>{0, eta});}
        void buildLinearKnots(const Box& box);
        void buildComplexKnots();
        void buildR0();
};






    
