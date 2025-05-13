#pragma once

#include "Input.h"
#include "Laser.h"


using lm_pair = std::pair<int,int>;
using lm_map = std::map<lm_pair,int>;
using block_map = std::map<int,lm_pair>;


class Angular
{
    public:




        explicit Angular(const Input& input) 
        : lmax{input.getJSON().at("Angular").at("lmax")}
        , mmin{input.getJSON().at("Angular").at("mmin")}
        , mmax{input.getJSON().at("Angular").at("mmax")}
        {}

        int getLmax() const {return lmax;}
        int getMmin() const {return mmin;}
        int getMmax() const {return mmax;}
        int getNlm() const {return nlm;}
        const lm_map& getLMMap() const { return lm_to_block;}
        const block_map& getBlockMap() const{return block_to_lm;}

        void buildMaps(const Laser& laser, const std::array<int,3>& initial_state);
        void printConfiguration(int rank);
        void dumpTo(const std::string& directory,int rank);

        
        
        
        

    private:

        // Member List Initialized
        int lmax{};
        int mmin{};
        int mmax{};
        

        // Derived
        int nlm{};
        lm_map lm_to_block{};
        block_map block_to_lm{};

        // Member Functions
        void buildZ(int m_i);
        void buildXY(int l_i, int m_i);
        void buildXYZ();
        void buildOdd();
        void buildEven();
};






