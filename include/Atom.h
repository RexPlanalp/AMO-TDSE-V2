#pragma once

#include <string>
#include <complex>
#include <stdexcept>
#include "Input.h"

class Atom {
public:

    struct Potentials 
    {
        static std::complex<double> hydrogen(const std::complex<double>& x) 
        {
            return -std::complex<double>{1.0} / (x + std::complex<double>{1e-25});
        }

        static std::complex<double> hydrogenDeriv(const std::complex<double>& x) 
        {
            return  std::complex<double>{1.0} / (x * x + std::complex<double>{1e-25});
        }
    };



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

    void printConfiguration(int rank);

private:
    std::complex<double> (*potential)(const std::complex<double>&);
    std::complex<double> (*derivative)(const std::complex<double>&);
    std::string    species;
};
