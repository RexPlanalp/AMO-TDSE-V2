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
        Basis() = default;


        explicit Basis(const Input& input,const Box& box)
        : m_Nbasis{input.getJSON().at("Basis").at("nbasis")}
        , m_order{input.getJSON().at("Basis").at("order")}
        , m_R0{input.getJSON().at("Basis").at("R0r")}
        , m_eta{input.getJSON().at("Basis").at("etar")}
        , m_spacing{input.getJSON().at("Basis").at("spacing")}
        {
            buildParameters();
            buildKnots(box);
        }

        int getNbasis() const {return m_Nbasis;}
        int getOrder() const {return m_order;}
        int getDegree() const {return m_degree;}
        double getR0() const {return m_R0;}
        double getEta() const {return m_eta;}
        const std::string& getSpacing() const {return m_spacing;}
        const std::vector<std::complex<double>> getKnots() const {return m_knots;}
        const std::vector<std::complex<double>> getComplexKnots() const {return m_complexKnots;}
        const std::vector<double> getRoots() const {return m_roots;}
        const std::vector<double> getWeights() const {return m_weights;}
        
        std::complex<double> integrateMatrixElement(int i, int j,MatrixIntegrand integrand,bool use_ecs) const;

        


        
    

    private:
        // Member List Initialized
        int m_Nbasis{};
        int m_order{};
        double m_R0{};
        double m_eta{};
        std::string m_spacing{};

        // Derived
        int m_degree{};
        std::vector<std::complex<double>> m_knots{};
        std::vector<std::complex<double>> m_complexKnots{};
        std::vector<double> m_roots{};
        std::vector<double> m_weights{};

        // Private Methods
        std::complex<double> ecs_x(double x) const {return (x < getR0()) ? (std::complex<double>{x, 0.0}) : getR0() + (x - getR0()) * std::exp(std::complex<double>{0, getEta()});}
        std::complex<double> ecs_w(double x, double w) const { return (x < getR0()) ? (std::complex<double>{w, 0.0}) : w * std::exp(std::complex<double>{0, getEta()});}
        void buildParameters();
        void buildKnots(const Box& box);
        void buildLinearKnots(const Box& box);
        void buildComplexKnots();
        void buildR0();
};






    
