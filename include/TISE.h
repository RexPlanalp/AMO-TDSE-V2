#pragma once 


#include "PetscWrappers/PetscEPS.h"
#include "Input.h"
#include "Basis.h"
#include "Atom.h"
#include "Angular.h"

class TISE
{
    public:
        TISE() = default;

        explicit TISE(const Input& input) 
        : m_maxIter(input.getJSON().at("TISE").at("max_iter"))
        , m_tolerance(input.getJSON().at("TISE").at("tolerance"))
        , m_nmax(input.getJSON().at("TISE").at("nmax"))
        {}

        PetscInt getMaxIter() const {return m_maxIter;}
        PetscReal getTol() const {return m_tolerance;}
        PetscInt getNmax() const {return m_nmax;}

       
    
    private:
        // List Initialized
        PetscInt m_maxIter{};
        PetscReal m_tolerance{};
        PetscInt m_nmax{};
};


inline std::ostream& operator<<(std::ostream& out, const TISE& tise)
{
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    out << "TISE Configuration: " << "\n\n";
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
    out << "maxIter: " << tise.getMaxIter() << "\n\n";
    out << "Tolerance: " << tise.getTol() << "\n\n";
    out << "nmax: " << tise.getNmax() << "\n\n";

    return out;
}