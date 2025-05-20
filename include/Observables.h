#pragma once

#include "Input.h"

class Observables
{   
    public:
        Observables() = default;

        Observables(const Input& input)
        : m_projOutBound{input.getJSON().at("Observables").at("Block").at("projOutBound")}
        , m_Emin{input.getJSON().at("Observables").at("PES").at("Emin")}
        , m_Emax{input.getJSON().at("Observables").at("PES").at("Emax")}
        , m_slice{input.getJSON().at("Observables").at("PES").at("slice")}
        {
            buildNe();
        }

        // Build methods
        void buildNe();
        
        bool getProjOut() const {return m_projOutBound;}
        double getEmin() const {return m_Emin;}
        double getEmax() const {return m_Emax;}
        const std::string& getSlice() const {return m_slice;}
        int getNe() const {return m_Ne;}
 

  


    private:
        bool m_projOutBound{};
        double m_Emin{};
        double m_Emax{};
        int m_Ne{};
        std::string m_slice{};
};


inline std::ostream& operator<<(std::ostream& out, const Observables& observables) 
{ 
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    out << "Observables Configuration: " << "\n\n";
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    
    out << "projOutBound: " << observables.getProjOut() << "\n\n";
    out << "Emin: " << observables.getEmin() << "\n\n";
    out << "Emax: " << observables.getEmax() << "\n\n";
    out << "Ne: " << observables.getNe() << "\n\n";
    out << "Slice: " << observables.getSlice() << "\n\n";

    return out;
}