#pragma once

#include "nlohmann/json.hpp"
#include "misc.h"

struct Input
{   

    Input(const std::string& inputPath) 
    : params(loadJson(inputPath))
    {validate();}


    void validate() const {}
            
    nlohmann::json params;
};