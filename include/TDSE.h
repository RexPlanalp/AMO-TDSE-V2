#pragma once

#include "Input.h"
#include  "petsc.h"

class TDSE
{   
    public:
        TDSE() = default;


        explicit TDSE(const Input& input)
        : m_initialNLM{input.getJSON().at("TDSE").at("initialNLM").get<std::array<int,3>>()}
        , m_tolerance{input.getJSON().at("TDSE").at("tolerance")}
        , m_maxIter{input.getJSON().at("TDSE").at("maxIter")}
        , m_restart{input.getJSON().at("TDSE").at("restart")}
        {}

        int getInitialN() const {return m_initialNLM[0];}
        int getInitialL() const {return m_initialNLM[1];}
        int getInitialM() const {return m_initialNLM[2];}
        PetscReal getTol() const {return m_tolerance;}
        PetscInt getMaxIter() const {return m_maxIter;}
        PetscInt getRestart() const {return m_restart;}
        

       
    private:
        std::array<int,3> m_initialNLM{};
        PetscReal m_tolerance{};
        PetscInt m_maxIter{};
        PetscInt m_restart{};
};

inline std::ostream& operator<<(std::ostream& out, const TDSE& tdse)
{
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    out << "TDSE Configuration: " << "\n\n";
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    

    out << "maxIter: " << tdse.getMaxIter() << "\n\n";
    out << "Restart: " << tdse.getRestart() << "\n\n";
    out << "Tol: " << tdse.getTol() << "\n\n";
    out << "Initial nlm: " << "n = " << tdse.getInitialN()  << " l = " << tdse.getInitialL()  << " m = " << tdse.getInitialM() << "\n\n";

    return out;
}