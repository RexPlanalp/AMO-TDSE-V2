#pragma once 


#include "PetscWrappers/PetscEPS.h"
#include "Input.h"
#include "BSpline.h"
#include "Atom.h"
#include "Angular.h"

class TISE
{
    public:
        explicit TISE(const Input& input) 
        : maxIter(input.getJSON().at("TISE").at("max_iter"))
        , tolerance(input.getJSON().at("TISE").at("tolerance"))
        , status{input.getJSON().at("TISE").at("status")}
        , outputPath{input.getJSON().at("TISE").at("outputPath")}
        , nmax(input.getJSON().at("TISE").at("nmax"))
        {}

        PetscInt getMaxIter() const {return maxIter;}
        PetscReal getTol() const {return tolerance;}
        bool getStatus() const {return status;}
        const std::string& getOutputPath() const {return outputPath;}
        PetscInt getNmax() const {return nmax;}

        void solve(const BSpline& bspline, const Atom& atom, const Angular& angular);
        void printConfiguration(int rank);
    
    private:
        // Member List Initialized
        PetscInt maxIter{};
        PetscReal tolerance{};
        bool status{};
        std::string outputPath{};
        PetscInt nmax{};

        std::string eigenvalueGroup = "eigenvalues";
        std::string eigenvectorGroup = "eigenvectors";

       
        
    
    

};