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
