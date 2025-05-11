#pragma once

#include <string>
#include <complex>
#include <stdexcept>
#include "Input.h"
#include "Potentials.h"

class Atom {
public:

    



    explicit Atom(const Input& input)
      : species{input.getJSON().at("Atom").at("potential")}
    {
        if (species == "H") 
        {
            potential  = &Potentials::hydrogen;
            derivative = &Potentials::hydrogenDeriv;
        }
    }

    std::complex<double> operator()(const std::complex<double>& x) const {return (*potential)(x);}

    std::complex<double> dpotential(const std::complex<double>& x) const {return (*derivative)(x);}

    const std::string& getSpecies() const {return species;}

    void printConfiguration(int rank);

private:
    std::complex<double> (*potential)(const std::complex<double>&);
    std::complex<double> (*derivative)(const std::complex<double>&);
    std::string    species;
};
