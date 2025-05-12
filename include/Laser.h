#pragma once

#include "Input.h"
#include "misc.h"

class Laser
{   
    public:
        explicit Laser(const Input& input)
        : N_cycles{input.getJSON().at("Laser").at("N")}
        , dt{input.getJSON().at("Laser").at("timeSpacing")}
        , w{input.getJSON().at("Laser").at("w")}
        , I{input.getJSON().at("Laser").at("I")}
        , polarization{input.getJSON().at("Laser").at("polarization").get<std::array<double,3>>()}
        , poynting{input.getJSON().at("Laser").at("poynting").get<std::array<double,3>>()}
        , ell{input.getJSON().at("Laser").at("ell")}
        , cep{input.getJSON().at("Laser").at("cepr")}
        {
            I /= Constants::I_AU;
            A_0 = std::sqrt(I) / getW();

            t_max = getN() * 2 * M_PI / getW();
            cep *= M_PI;
            
            Nt = static_cast<int>(std::round(getTmax() / getTimeSpacing())) + 1;
            times.resize(Nt);
            for (int idx = 0; idx < Nt; ++idx)
            {
                times[idx] = idx * getTimeSpacing();
            }

            ellipticity = crossProduct(getPolarization(),getPoynting());

            normalize(polarization);
            normalize(poynting);
            normalize(ellipticity);
            buildNonzeroComponents();
        }

        double getN() const {return N_cycles;} 
        double getW() const {return w;} 
        double getTimeSpacing() const {return dt;} 
        double getI() const {return I;} 
        double getA0() const {return A_0;} 
        double getTmax() const {return t_max;}
        double getEll() const {return ell;}
        double getCEP() const {return cep;}
        int getNt() const {return Nt;}
        const std::array<double,3>& getPolarization() const {return polarization;}
        const std::array<double,3>& getPoynting() const {return poynting;}
        const std::array<double,3>& getEllipticity() const {return ellipticity;}
        const std::array<int,3>& getComponents() const {return components;}

        void printConfiguration(int rank);
        double getTime(int i) const;
        void dumpTo(const std::string& directory,int rank);
        double A(double t, int idx) const;




    
    private:  
        // Member List initialized
        double N_cycles{};
        double dt{};
        double w{};
        double I{};
        std::array<double,3> polarization{};
        std::array<double,3> poynting{};
        double ell{};
        double cep{};

        // Derived
        double A_0{};
        double t_max{};
        int Nt{};
        std::array<double,3> ellipticity{};
        std::array<int,3> components{};
        std::vector<double> times{};

        // Member functions
        void buildNonzeroComponents();
        double sin2_envelope(double t) const;
        void validateInput();

};