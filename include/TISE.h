#pragma once 

#include "nlohmann/json.hpp"
#include "PetscWrappers/PetscEPS.h"

class TISE
{
    public:
        TISE() = delete;

        explicit TISE(const nlohmann::json& input_file) 
        : maxIter(input_file.at("TISE").at("max_iter"))
        , tolerance(input_file.at("TISE").at("tolerance"))
        , on{input_file.at("TISE").at("on")}
        , outputPath{input_file.at("TISE").at("output")}
        , nmax(input_file.at("TISE").at("n_max"))
        {}
    
    private:
        // Member List Initialized
        PetscInt maxIter{};
        PetscReal tolerance{};
        bool on{};
        std::string outputPath{};
        PetscInt nmax{};

        
    
    

};