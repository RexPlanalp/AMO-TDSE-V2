#pragma once

#include <string>
#include <complex>
#include <stdexcept>
#include "Input.h"
#include "Potentials.h"
#include "MatrixElements.h"

class Atom {
    public:

        Atom() = default;

        explicit Atom(const Input& input)
        : species{input.getJSON().at("Atom").at("potential")}
        {
            if (species == "H") 
            {
                potential  = &Potentials::hydrogen;
                derivative = &Potentials::hydrogenDeriv;
                potentialType = RadialMatrixType::H;
            }
        }

        std::complex<double> operator()(const std::complex<double>& x) const {return (*potential)(x);}

        std::complex<double> dpotential(const std::complex<double>& x) const {return (*derivative)(x);}

        const std::string& getSpecies() const {return species;}

        RadialMatrixType getType() const {return potentialType;}

        void printConfiguration(int rank);

    private:
        std::complex<double> (*potential)(const std::complex<double>&);
        std::complex<double> (*derivative)(const std::complex<double>&);
        std::string    species;
        RadialMatrixType potentialType;
};
