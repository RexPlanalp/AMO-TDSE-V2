#pragma once

#include "Input.h"

class Box 
{
    public:

        explicit Box(const Input& input) 
        : gridSize{input.getJSON().at("Box").at("gridSize")}
        , gridSpacing{input.getJSON().at("Box").at("gridSpacing")}
        {   
           

            Nr = static_cast<int>(std::round(getGridSize() / getGridSpacing())) + 1;

            positions.resize(Nr);
            for (int idx = 0; idx < Nr; ++idx)
            {
                positions[idx] = idx * getGridSpacing();
            }

        

        }

        double getGridSize() const {return gridSize;}
        double getGridSpacing() const {return gridSpacing;}
        int getNr() const {return Nr;}

        double getPosition(int i) const;
        void printConfiguration(int rank) const;

    private:
        // Member List Initialized
        double gridSize{};
        double gridSpacing{};
        
        // Derived
        int Nr{};
        std::vector<double> positions{};

        // Member Functions
        void validateInput();
};