#pragma once

#include "misc.h"

class Input
{   
    public:
        Input(const std::string& inputPath) 
        : params(loadJson(inputPath))
        {}

        const nlohmann::json& getJSON() const {return params;}


    private:     
        nlohmann::json params;
};