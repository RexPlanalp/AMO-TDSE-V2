#pragma once

#include <string>
#include "nlohmann/json.hpp"
#include <optional>
#include "Laser.h"

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

        void buildMaps(const Laser& laser, const std::array<int,3>& initial_state);
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



namespace AngularElement
{
    inline double f(int l, int m)
    {   
        int numerator = (l+1)*(l+1) - m*m;
        int denominator = (2*l + 1)*(2*l+3);
        return sqrt(numerator/double(denominator));
    }

    inline double g(int l, int m)
    {
        int numerator = l*l - m*m;
        int denominator = (2*l-1)*(2*l+1);
        return sqrt(numerator/double(denominator));
    }

    inline double a(int l, int m)
    {
        int numerator = (l+m);
        int denominator = (2*l +1) * (2*l-1);
        double f1 = sqrt(numerator/double(denominator));
        double f2 = - m * std::sqrt(l+m-1) - std::sqrt((l-m)*(l*(l-1)-m*(m-1)));
        return f1*f2;
        
    }

    inline double atilde(int l, int m)
    {
        int numerator = (l-m);
        int denominator = (2*l+1)*(2*l-1);
        double f1 = sqrt(numerator/double(denominator));
        double f2 = - m * std::sqrt(l-m-1) + std::sqrt((l+m)*(l*(l-1)-m*(m+1)));
        return f1*f2;
    }

    inline double b(int l, int m)
    {
        return -atilde(l+1,m-1);
    }

    inline double btilde(int l, int m)
    {
        return -a(l+1,m+1);
    }

    inline double d(int l, int m)
    {
        double numerator = (l-m+1)*(l-m+2);
        double denominator = (2*l+1)*(2*l+3);
        return std::sqrt(numerator/double(denominator));
    }

    inline double dtilde(int l, int m)
    {
        return d(l,-m);
    }

    inline double c(int l, int m)
    {
        return dtilde(l-1,m-1);
    }

    inline double ctilde(int l, int m)
    {
        return d(l-1,m+1);
    }

    inline double alpha(int l, int m)
    {
        double numerator = (l+m-1)*(l+m);
        double denominator = 4*(2*l+1)*(2*l-1);
        return std::sqrt(numerator/denominator);
    }

    inline double beta(int l, int m)
    {
        double numerator = (l-m+1)*(l-m+2)*(l+1);
        double denominator = 2*(2*l+1)*(2*l+2)*(2*l+3);
        return -std::sqrt(numerator/denominator);
    }

    inline double charlie(int l, int m)
    {
        double numerator = (l-m-1)*(l-m);
        double denominator = 4*(2*l+1)*(2*l-1);
        return std::sqrt(numerator/denominator);
    }

    inline double delta(int l, int m)
    {
        double numerator = (l+m+1)*(l+m+2)*(l+1);
        double denominator = 2*(2*l+1)*(2*l+2)*(2*l+3);
        return -std::sqrt(numerator/denominator);
    }

    inline double echo(int l, int m)
    {
        double numerator = (l+m)*(l-m);
        double denominator = (2*l-1)*(2*l+1);
        return std::sqrt(numerator/denominator);
    }

    inline double foxtrot(int l, int m)
    {
        double numerator = (l+m+1)*(l-m+1);
        double denominator = (2*l+1)*(2*l+3);
        return std::sqrt(numerator/denominator);
    }
}


