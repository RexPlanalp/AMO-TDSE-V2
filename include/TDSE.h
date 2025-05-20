#pragma once

#include "Input.h"
#include "petsc.h"
#include "PetscWrappers/PetscVec.h"
#include "TISE.h"
#include "Basis.h"


#include "PetscWrappers/PetscMat.h"

class TDSE
{   
    public:
        TDSE() = default;


        explicit TDSE(const Input& input)
        : m_status{input.getJSON().at("TDSE").at("status")}
        , m_initialNLM{input.getJSON().at("TDSE").at("initialNLM").get<std::array<int,3>>()}
        , m_tolerance{input.getJSON().at("TDSE").at("tolerance")}
        , m_maxIter{input.getJSON().at("TDSE").at("maxIter")}
        , m_restart{input.getJSON().at("TDSE").at("restart")}
        {}

        bool getStatus() const {return m_status;}
        int getInitialN() const {return m_initialNLM[0];}
        int getInitialL() const {return m_initialNLM[1];}
        int getInitialM() const {return m_initialNLM[2];}
        PetscReal getTol() const {return m_tolerance;}
        PetscInt getMaxIter() const {return m_maxIter;}
        PetscInt getRestart() const {return m_restart;}
        

       
    private:
        bool m_status{};
        std::array<int,3> m_initialNLM{};
        PetscReal m_tolerance{};
        PetscInt m_maxIter{};
        PetscInt m_restart{};
};