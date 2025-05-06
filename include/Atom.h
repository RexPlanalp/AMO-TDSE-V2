#include <complex>
#include "Potentials.h"
#include "nlohmann/json.hpp"

class Atom
{   public:
        Atom() = delete;
        explicit Atom(const nlohmann::json& input_file)
        : species(input_file.at("Atom").at("potential"))
        {}

        template<typename T>
        T potentialVal(T position)
        {
        if (species == "H")
            return hydrogenPotential(position);

        throw std::runtime_error("unknown species");
        }

        template<typename T>
        T potentialDerivativeVal(T position)
        {
        if (species == "H")
            return hydrogenPotentialDerivative(position);

        throw std::runtime_error("unknown species");
        }

    private:
        std::string species{};
};


