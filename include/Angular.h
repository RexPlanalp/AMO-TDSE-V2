#pragma once

#include <string>
#include "nlohmann/json.hpp"
#include <optional>

using lm_pair = std::pair<int,int>;
using lm_map = std::map<lm_pair,int>;
using block_map = std::map<int,lm_pair>;


class Angular
{
    public:

        Angular() = delete;

        explicit Angular(const nlohmann::json& input_file) 
        : l_max{input_file.at("Angular").at("l_max")}
        , m_min{input_file.at("Angular").at("m_min")}
        , m_max{input_file.at("Angular").at("m_max")}
        {
            if (LMax() <= 0)
            {
                throw std::invalid_argument("l_max size must be greater than zero. You entered: " + std::to_string(LMax()));
            }
            

            if ((abs(MMin()) > LMax()) || (abs(MMax()) > LMax() ))
            {
                throw std::invalid_argument("m_min and m_max should be less than or equal in magnitude to l_max.");
            }
            
        }

        int LMax() const {return l_max;}
        int MMin() const {return m_min;}
        int MMax() const {return m_max;}
        int N_lm() const {return n_lm;}

        const lm_map& LMMap()  
        {
            if (!lm_to_block.has_value())
            {
                throw std::logic_error("lm_to_block has not been built yet.");
            }
            

            return lm_to_block.value();
        }

        const block_map& BlockMap() const
        {
            if (!block_to_lm.has_value())
            {
                throw std::logic_error("block_to_lm has not been built yet.");
            }
            return block_to_lm.value();
        }

        void buildMaps(const std::array<int,3>& components, const std::array<int,3>& initial_state);
        void buildZ(int m_i);
        void buildXY(int l_i, int m_i);
        void buildXYZ();
        void buildOdd();
        void buildEven();
        void dumpTo(const std::string& directory,int rank);

        
        
        
        

    private:
        int l_max{};
        int m_min{};
        int m_max{};
        int n_lm{};

        std::optional<lm_map> lm_to_block{};
        std::optional<block_map> block_to_lm{};


};