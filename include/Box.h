#pragma once

#include <string>
#include "nlohmann/json.hpp"

class Box 
{
    public:

        Box() = delete;

        explicit Box(const nlohmann::json& input_file) 
        : grid_size{input_file["Box"]["grid_size"].get<double>()}
        , grid_spacing{input_file["Box"]["grid_spacing"].get<double>()}
        {validateInput();}

        double GridSize() const {return grid_size;}
        double GridSpacing() const {return grid_spacing;}

        int Nr() const;
        double Position(int i) const;

    private:
        // Member List Initialized
        double grid_size{};
        double grid_spacing{};

        // Member Functions
        void validateInput();
};