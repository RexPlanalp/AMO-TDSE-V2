#pragma once

#include "Input.h"
#include "misc.h"

class Laser
{   
    public:
        Laser() = default;


        explicit Laser(const Input& input)
        : m_Ncycles{input.getJSON().at("Laser").at("N")}
        , m_dt{input.getJSON().at("Laser").at("timeSpacing")}
        , m_w{input.getJSON().at("Laser").at("w")}
        , m_I{input.getJSON().at("Laser").at("I")}
        , m_polarization{input.getJSON().at("Laser").at("polarization").get<std::array<double,3>>()}
        , m_poynting{input.getJSON().at("Laser").at("poynting").get<std::array<double,3>>()}
        , m_ell{input.getJSON().at("Laser").at("ell")}
        , m_cep{input.getJSON().at("Laser").at("cepr")}
        {
            buildParameters();
            buildNt();
            buildTimes();
            buildVectors();
        }

        void dumpTo(const std::string& directory);

        // Getters
        double getN() const {return m_Ncycles;} 
        double getW() const {return m_w;} 
        double getTimeSpacing() const {return m_dt;} 
        double getI() const {return m_I;} 
        double getA0() const {return m_A_0;} 
        double getTmax() const {return m_Tmax;}
        double getEll() const {return m_ell;}
        double getCEP() const {return m_cep;}
        int getNt() const {return m_Nt;}
        const std::array<double,3>& getPolarization() const {return m_polarization;}
        const std::array<double,3>& getPoynting() const {return m_poynting;}
        const std::array<double,3>& getEllipticity() const {return m_ellipticity;}
        const std::array<int,3>& getComponents() const {return m_components;}
        double operator[](int i) const;


        // Public Methods
        double A(double t, int idx) const;

        



    
    private:  
        // Member List initialized
        double m_Ncycles{};
        double m_dt{};
        double m_w{};
        double m_I{};
        std::array<double,3> m_polarization{};
        std::array<double,3> m_poynting{};
        double m_ell{};
        double m_cep{};

        // Derived
        double m_A_0{};
        double m_Tmax{};
        int m_Nt{};
        std::array<double,3> m_ellipticity{};
        std::array<int,3> m_components{};
        std::vector<double> m_times{};

        //  Private  Member functions
        double sin2_envelope(double t) const;
        void validateInput();
        void buildParameters();
        void buildNt();
        void buildTimes();
        void buildVectors();

};



inline std::ostream& operator<<(std::ostream& out, const Laser& laser)
{

        out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        out << "Laser Configuration: " << "\n\n";
        out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
        out << "Ncycles: " << laser.getN() << "\n\n";
        out << "dt: " << laser.getTimeSpacing() << "\n\n";
        out << "tmax: " << laser.getTmax() << "\n\n";
        out << "Nt: " << laser.getNt() << "\n\n";
        out << "w: " << laser.getW() << "\n\n";
        out << "I: " << laser.getI() << "\n\n";
        out << "A_0: " << laser.getA0() << "\n\n";
        out << "ell: " << laser.getEll() << "\n\n";
        out << "CEP: " << laser.getCEP() << "\n\n";
        out << "polarization: " << laser.getPolarization()[0] << " " << laser.getPolarization()[1] << " " << laser.getPolarization()[2] <<  "\n\n";
        out << "poynting: " << laser.getPoynting()[0] << " " << laser.getPoynting()[1] << " " << laser.getPoynting()[2] <<  "\n\n";
        out << "ellipticity: " << laser.getEllipticity()[0] << " " << laser.getEllipticity()[1] << " " << laser.getEllipticity()[2] <<  "\n\n";
        out << "components: " << laser.getComponents()[0] << " " << laser.getComponents()[1] << " " << laser.getComponents()[2] <<  "\n\n";

        return out;
}