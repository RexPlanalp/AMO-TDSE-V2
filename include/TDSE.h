#pragma once

#include "Input.h"
#include "petsc.h"
#include "PetscWrappers/PetscVec.h"
#include "TISE.h"
#include "Basis.h"
#include "Angular.h"

#include "PetscWrappers/PetscMat.h"

class TDSE
{   public:
        explicit TDSE(int rank,const Input& input)
        : status{input.getJSON().at("TDSE").at("status")}
        , outputPath{input.getJSON().at("TDSE").at("outputPath")}
        , initialNLM{input.getJSON().at("TDSE").at("initialNLM").get<std::array<int,3>>()}
        , tolerance{input.getJSON().at("TDSE").at("tolerance")}
        , maxIter{input.getJSON().at("TDSE").at("maxIter")}
        , restart{input.getJSON().at("TDSE").at("restart")}
        , hhgStatus{input.getJSON().at("TDSE").at("HHG")}
        {
            createDirectory("TDSE",rank);
        }

        bool getStatus() const {return status;}
        const std::string& getOutputPath() const {return outputPath;}
        const std::array<int,3>& getInitialNLM() const {return initialNLM;}
        PetscReal getTol() const {return tolerance;}
        PetscInt getMaxIter() const {return maxIter;}
        PetscInt getRestart() const {return restart;}
        bool getHHGStatus() const {return hhgStatus;}

        void solve(int rank,const TISE& tise,const Basis& basis, const Angular& angular, const Atom& atom, const Laser& laser);
        void printConfiguration(int rank);

        Vector loadInitialState(const TISE& tise,const Basis& basis, const Angular& angular);
        std::pair<Matrix,Matrix> constructAtomicInteraction(const Basis& basis, const Angular& angular,const Atom& atom, const Laser& laser);
        
        Matrix constructZInteraction(const Basis& basis, const Angular& angular);
        std::pair<Matrix,Matrix> constructXYInteraction(const Basis& basis, const Angular& angular);

        Matrix constructXHHG(const Basis& basis, const Angular& angular);
        Matrix constructYHHG(const Basis& basis, const Angular& angular);
        Matrix constructZHHG(const Basis& basis, const Angular& angular);

        Matrix constructAtomicS(const Basis& basis, const Angular& angular);

    private:
        Matrix kroneckerProduct(const Matrix& A, const Matrix& B, PetscInt nnz_A, PetscInt nnz_B);
    
        bool status{};
        std::string outputPath{};
        std::array<int,3> initialNLM{};
        PetscReal tolerance{};
        PetscInt maxIter{};
        PetscInt restart{};
        bool hhgStatus{};

        std::string outputGroup = "";
        std::string outputName = "psiFinal";
};