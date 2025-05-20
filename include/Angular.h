#pragma once

#include "Input.h"
#include "Laser.h"
#include "TDSE.h"

using lmPair = std::pair<int,int>;
using lmMap = std::map<lmPair,int>;
using blockMap = std::map<int,lmPair>;


class Angular
{
    public:

        Angular() = default;


        explicit Angular(const Input& input, const Laser& laser, const TDSE& tdse) 
        : m_lmax{input.getJSON().at("Angular").at("lmax")}
        , m_mmin{input.getJSON().at("Angular").at("mmin")}
        , m_mmax{input.getJSON().at("Angular").at("mmax")}
        {
            buildMaps(laser,tdse);
        }

        

        // Getters
        int getLmax() const {return m_lmax;}
        int getMmin() const {return m_mmin;}
        int getMmax() const {return m_mmax;}
        int getNblocks() const {return m_Nblocks;}
        const lmMap& getLMMap() const {return m_lmToBlock;}
        const blockMap& getBlockMap() const{return m_blockTolm;}



        
        
        
        

    private:

        // List Initialized
        int m_lmax{};
        int m_mmin{};
        int m_mmax{};
        
        // Derived
        int m_Nblocks{};
        lmMap m_lmToBlock{};
        blockMap m_blockTolm{};

        // Private Member Functions
        void buildZ(int m_i);
        void buildXY(int l_i, int m_i);
        void buildXYZ();
        void buildOdd();
        void buildEven();
        void buildMaps(const Laser& laser, const TDSE& tdse);
};



inline std::ostream& operator<<(std::ostream& out, const Angular& angular)
{
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    out << "Angular Configuration: " << "\n\n";
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";

    out << "lmax: " << angular.getLmax() << "\n\n";
    out << "mmax: " << angular.getMmax() << "\n\n";
    out << "mmin: " << angular.getMmin() << "\n\n";
    out << "nlm: "  << angular.getNblocks() << "\n\n";
}





