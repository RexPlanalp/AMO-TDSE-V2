#pragma once

#include "Input.h"

class Box 
{
    public:
        Box() = default;


        explicit Box(const Input& input) 
        : m_gridSize{input.getJSON().at("Box").at("gridSize")}
        , m_gridSpacing{input.getJSON().at("Box").at("gridSpacing")}
        {   
            buildNr();
            buildPositions();
        }



        // Getters
        double getGridSize() const {return m_gridSize;}
        double getGridSpacing() const {return m_gridSpacing;}
        int getNr() const {return m_Nr;}
        double operator[](int i) const;

    private:
        // List Initialized
        double m_gridSize{};
        double m_gridSpacing{};

        // Derived
        int m_Nr{};
        std::vector<double> m_positions{};

        // Private Methods
        void buildNr();
        void buildPositions();
};

inline std::ostream& operator<<(std::ostream& out, const Box& box) 
{ 
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    out << "Box Configuration: " << "\n\n";
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    
    out << "rmax: " << box.getGridSize() << "\n\n";
    out << "dr: " << box.getGridSpacing() << "\n\n";
    out << "Nr: " << box.getNr() << "\n\n";
}