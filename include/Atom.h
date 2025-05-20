#pragma once


#include "Input.h"
#include "Potentials.h"
#include "MatrixElements.h"

class Atom {
    public:

        Atom() = default;

        explicit Atom(const Input& input)
        : m_species{input.getJSON().at("Atom").at("species")}
        {
            buildPotential();
        }

        // Getters
        std::complex<double> potential(const std::complex<double>& x) const {return (*m_potential)(x);}
        std::complex<double> derivative(const std::complex<double>& x) const {return (*m_derivative)(x);}
        const std::string& getSpecies() const {return m_species;}

    private:
        // List Initialized
        std::string m_species;

        // Derived
        std::complex<double> (*m_potential)(const std::complex<double>&);
        std::complex<double> (*m_derivative)(const std::complex<double>&);

        // Private Methods
        void buildPotential();
};
