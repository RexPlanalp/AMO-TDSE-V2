#pragma once

#include "nlohmann/json.hpp"
#include "misc.h"

class Input
{   
    public:
        Input(int rank,const std::string& inputPath) 
        : params(loadJson(inputPath))
        {
            createDirectory("misc",rank);
            createDirectory("images",rank);
        }

        const nlohmann::json& getJSON() const {return params;}


    private:     
        nlohmann::json params;
};