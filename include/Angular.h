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

    public:
        struct Coupling
        {
            static double f(int l, int m)
            {   
                int numerator = (l+1)*(l+1) - m*m;
                int denominator = (2*l + 1)*(2*l+3);
                return sqrt(numerator/double(denominator));
            }

            static double g(int l, int m)
            {
                int numerator = l*l - m*m;
                int denominator = (2*l-1)*(2*l+1);
                return sqrt(numerator/double(denominator));
            }

            static double a(int l, int m)
            {
                int numerator = (l+m);
                int denominator = (2*l +1) * (2*l-1);
                double f1 = sqrt(numerator/double(denominator));
                double f2 = - m * std::sqrt(l+m-1) - std::sqrt((l-m)*(l*(l-1)-m*(m-1)));
                return f1*f2;
                
            }

            static double atilde(int l, int m)
            {
                int numerator = (l-m);
                int denominator = (2*l+1)*(2*l-1);
                double f1 = sqrt(numerator/double(denominator));
                double f2 = - m * std::sqrt(l-m-1) + std::sqrt((l+m)*(l*(l-1)-m*(m+1)));
                return f1*f2;
            }

            static double b(int l, int m)
            {
                return -atilde(l+1,m-1);
            }

            static double btilde(int l, int m)
            {
                return -a(l+1,m+1);
            }

            static double d(int l, int m)
            {
                double numerator = (l-m+1)*(l-m+2);
                double denominator = (2*l+1)*(2*l+3);
                return std::sqrt(numerator/double(denominator));
            }

            static double dtilde(int l, int m)
            {
                return d(l,-m);
            }

            static double c(int l, int m)
            {
                return dtilde(l-1,m-1);
            }

            static double ctilde(int l, int m)
            {
                return d(l-1,m+1);
            }

            static double alpha(int l, int m)
            {
                double numerator = (l+m-1)*(l+m);
                double denominator = 4*(2*l+1)*(2*l-1);
                return std::sqrt(numerator/denominator);
            }

            static double beta(int l, int m)
            {
                double numerator = (l-m+1)*(l-m+2)*(l+1);
                double denominator = 2*(2*l+1)*(2*l+2)*(2*l+3);
                return -std::sqrt(numerator/denominator);
            }

            static double charlie(int l, int m)
            {
                double numerator = (l-m-1)*(l-m);
                double denominator = 4*(2*l+1)*(2*l-1);
                return std::sqrt(numerator/denominator);
            }

            static double delta(int l, int m)
            {
                double numerator = (l+m+1)*(l+m+2)*(l+1);
                double denominator = 2*(2*l+1)*(2*l+2)*(2*l+3);
                return -std::sqrt(numerator/denominator);
            }

            static double echo(int l, int m)
            {
                double numerator = (l+m)*(l-m);
                double denominator = (2*l-1)*(2*l+1);
                return std::sqrt(numerator/denominator);
            }

            static double foxtrot(int l, int m)
            {
                double numerator = (l+m+1)*(l-m+1);
                double denominator = (2*l+1)*(2*l+3);
                return std::sqrt(numerator/denominator);
            }
        };


};






