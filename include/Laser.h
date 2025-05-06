#pragma once

#include "nlohmann/json.hpp"
#include "misc.h"

class Laser
{   
    public:
        Laser() = delete;

        explicit Laser(const nlohmann::json& input_file)
        : N_cycles{input_file.at("Laser").at("N")}
        , w{input_file.at("Laser").at("w")}
        , dt{input_file.at("Laser").at("dt")}
        , I_si{input_file.at("Laser").at("I")}
        , polarization{input_file.at("Laser").at("polarization").get<std::array<double,3>>()}
        , poynting{input_file.at("Laser").at("poynting").get<std::array<double,3>>()}
        , ell{input_file.at("Laser").at("ell")}
        , cep{input_file.at("Laser").at("cep_r")}
        {
            I_au = I_si / 3.51E16;
            A_0 = std::sqrt(I_au) / W();

            t_max = N() * 2 * M_PI / W();
            cep *= M_PI;

            ellipticity = crossProduct(Polarization(),Poynting());

            normalize(polarization);
            normalize(poynting);
            normalize(ellipticity);

            nonzeroComponents();



            if (N() <= 0.0)
            {
                throw std::invalid_argument("N_cycles must be greater than or equal to zero. You entered: " + std::to_string(N()));
            }
            if (W() <= 0.0)
            {
                throw std::invalid_argument("w must be greater than or equal to zero. You entered: " + std::to_string(W()));
            }
            if (TimeSpacing() <= 0.0)
            {
                throw std::invalid_argument("dt must be greater than or equal to zero. You entered: " + std::to_string(TimeSpacing()));
            }
            if (realDotProduct(Polarization(),Poynting()) != 0.0)
            {
                throw std::invalid_argument("Polarization and Poynting vectors must be orthogonal.");
            }
            if (!((ELL() <= 1.0) && (ELL() >= 0.0)))
            {
                throw std::invalid_argument("Ell must be between 0.0 and 1.0. You entered: " + std::to_string(ELL()));
            }


        }

        double N() const {return N_cycles;} 
        double W() const {return w;} 
        double TimeSpacing() const {return dt;} 
        double I() const {return I_au;} 
        double A0() const {return A_0;} 
        double TMAX() const {return t_max;}

        const std::array<double,3>& Polarization() const {return polarization;}
        const std::array<double,3>& Poynting() const {return poynting;}
        const std::array<double,3>& Ellipticity() const {return ellipticity;}

        double ELL() const {return ell;}
        double CEP() const {return cep;}

        double sin2_envelope(double t);
        double A(double t, int idx);

        void dumpTo(const std::string& directory,int rank);

        int Nt() const;
        double Time(int i) const;



    
    private:  
        void nonzeroComponents();
        
        double N_cycles{};
        double w{};
        double dt{};
        double I_au{};
        double I_si{};
        double A_0{};
        double t_max{};

        std::array<double,3> polarization{};
        std::array<double,3> poynting{};
        std::array<double,3> ellipticity{};
        std::array<int,3> components{};
        double ell{};
        double cep{};

};