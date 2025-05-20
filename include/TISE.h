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
        , m_status{input.getJSON().at("TISE").at("status")}
        , m_nmax(input.getJSON().at("TISE").at("nmax"))
        {}

        PetscInt getMaxIter() const {return m_maxIter;}
        PetscReal getTol() const {return m_tolerance;}
        bool getStatus() const {return m_status;}
        PetscInt getNmax() const {return m_nmax;}

       
    
    private:
        // List Initialized
        PetscInt m_maxIter{};
        PetscReal m_tolerance{};
        bool m_status{};
        PetscInt m_nmax{};
};