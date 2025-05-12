#pragma once

#include "Input.h"
#include "petsc.h"

#include "PetscWrappers/PetscMat.h"

class TDSE
{   public:
        explicit TDSE(const Input& input)
        : status{input.getJSON().at("TDSE").at("status")}
        , outputPath{input.getJSON().at("TDSE").at("outputPath")}
        , initialState{input.getJSON().at("TDSE").at("initialState").get<std::array<int,3>>()}
        , tolerance{input.getJSON().at("TDSE").at("tolerance")}
        , maxIter{input.getJSON().at("TDSE").at("maxIter")}
        , krylovDim{input.getJSON().at("TDSE").at("krylovDim")}
        {}

        bool getStatus() const {return status;}
        const std::string& getOutputPath() const {return outputPath;}
        const std::array<int,3>& getInitialState() const {return initialState;}
        PetscReal getTol() const {return tolerance;}
        PetscInt getMaxIter() const {return maxIter;}
        PetscInt getKrylovDim() const {return krylovDim;}

        void printConfiguration(int rank);

    private:
        Matrix kroneckerProduct(const Matrix& A, const Matrix& B, PetscInt nnz_A, PetscInt nnz_B);
    
        bool status{};
        std::string outputPath{};
        std::array<int,3> initialState{};
        PetscReal tolerance{};
        PetscInt maxIter{};
        PetscInt krylovDim{};
};