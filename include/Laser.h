#pragma once

#include "nlohmann/json.hpp"
#include "misc.h"

class Laser
{   
    public:
        Laser() = delete;

        explicit Laser(const nlohmann::json& input_file)
        : N_cycles{input_file.at("Laser").at("N")}
        , dt{input_file.at("Laser").at("time_spacing")}
        , w{input_file.at("Laser").at("w")}
        , I_si{input_file.at("Laser").at("I")}
        , polarization{input_file.at("Laser").at("polarization").get<std::array<double,3>>()}
        , poynting{input_file.at("Laser").at("poynting").get<std::array<double,3>>()}
        , ell{input_file.at("Laser").at("ell")}
        , cep_r{input_file.at("Laser").at("cep_r")}
        {
            validateInput();

            I_au = I_si / 3.51E16;
            A_0 = std::sqrt(I_au) / W();

            t_max = N() * 2 * M_PI / W();
            cep = cep_r * M_PI;

            ellipticity = crossProduct(Polarization(),Poynting());

            normalize(polarization);
            normalize(poynting);
            normalize(ellipticity);
            buildNonzeroComponents();
        }

        double N() const {return N_cycles;} 
        double W() const {return w;} 
        double TimeSpacing() const {return dt;} 
        double I() const {return I_au;} 
        double A0() const {return A_0;} 
        double TMAX() const {return t_max;}
        double ELL() const {return ell;}
        double CEP() const {return cep;}
        double CEP_R() const {return cep_r;}
        const std::array<double,3>& Polarization() const {return polarization;}
        const std::array<double,3>& Poynting() const {return poynting;}
        const std::array<double,3>& Ellipticity() const {return ellipticity;}
        const std::array<int,3>& Components() const {return components;}

        int Nt() const;
        double Time(int i) const;
        void dumpTo(const std::string& directory,int rank);




    
    private:  
        // Member List initialized
        double N_cycles{};
        double dt{};
        double w{};
        double I_si{};
        std::array<double,3> polarization{};
        std::array<double,3> poynting{};
        double ell{};
        double cep_r{};

        // Default Initialized
        double A_0;
        double I_au;
        double t_max;
        std::array<double,3> ellipticity;
        std::array<int,3> components;
        double cep;

        // Member functions
        void buildNonzeroComponents();
        double sin2_envelope(double t);
        double A(double t, int idx);
        void validateInput();

};