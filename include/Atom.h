#pragma once

#include <complex>
#include "Potentials.h"
#include "nlohmann/json.hpp"

class Atom
{   public:
        Atom() = delete;
        explicit Atom(const nlohmann::json& input_file)
        : species(input_file.at("Atom").at("potential"))
        {validateInput();}

        template<typename T>
        T potentialVal(T position)
        {
        if (Species() == "H")
            return Potentials::hydrogenPotential(position);

        throw std::runtime_error("unknown species");
        }

        template<typename T>
        T potentialDerivativeVal(T position)
        {
        if (Species() == "H")
            return Potentials::hydrogenPotentialDerivative(position);

        throw std::runtime_error("unknown species");
        }

        const std::string& Species() const {return species;}

    private:

        // Member List Initialized
        std::string species{};

        // Member Functions
        void validateInput();
};


