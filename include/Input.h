#pragma once

#include "nlohmann/json.hpp"
#include "misc.h"

class Input
{   
    public:
        Input(const std::string& inputPath) 
        : params(loadJson(inputPath))
        {}


        void validate() const {}
        const nlohmann::json& getJSON() const {return params;}


    private:     
        nlohmann::json params;
};