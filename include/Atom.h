#pragma once


#include "Input.h"
#include "Potentials.h"
#include "MatrixElements.h"

enum class RadialMatrixType;

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
        RadialMatrixType getType() const {return m_type;}

    private:
        // List Initialized
        std::string m_species{};

        // Derived
        std::complex<double> (*m_potential)(const std::complex<double>&);
        std::complex<double> (*m_derivative)(const std::complex<double>&);
        RadialMatrixType m_type{};

        // Private Methods
        void buildPotential();
};


inline std::ostream& operator<<(std::ostream& out,const Atom& atom)
{
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    out << "Atom Configuration: " << "\n\n";
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";

    out << "Species: " << atom.getSpecies() <<  "\n\n";

    return out;
}