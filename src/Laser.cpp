#include "Laser.h"

void Laser::validateInput()
{
    if (N() <= 0.0)
    {
        throw std::invalid_argument("N_cycles must be greater than or equal to zero. You entered: " + std::to_string(N()));
    }
    if (TimeSpacing() <= 0.0)
    {
        throw std::invalid_argument("dt must be greater than or equal to zero. You entered: " + std::to_string(TimeSpacing()));
    }
    if (W() <= 0.0)
    {
        throw std::invalid_argument("w must be greater than or equal to zero. You entered: " + std::to_string(W()));
    }
    if (I() <= 0.0)
    {
        throw std::invalid_argument("I must be greater than or equal to zero. You entered: " + std::to_string(I()));
    }
    if (!Polarization()[0] && !Polarization()[1] && !Polarization()[2])
    {
        throw std::invalid_argument("Polarization must have a non-zero component.");
    }
    if (!Poynting()[0] && !Poynting()[1] && !Poynting()[2])
    {
        throw std::invalid_argument("Polarization must have a non-zero component.");
    }
    if (realDotProduct(Polarization(),Poynting()) != 0.0)
    {
        throw std::invalid_argument("Polarization and Poynting vectors must be orthogonal.");
    }
    if (!((CEP_R() >= 0.0) && (CEP_R() <= 1.0)))
    {
        throw std::invalid_argument("CEP ratio should be between 0.0 and 1.0. You entered: " + std::to_string(CEP_R()));
    }
    if (!((ELL() <= 1.0) && (ELL() >= 0.0)))
    {
        throw std::invalid_argument("Ell must be between 0.0 and 1.0. You entered: " + std::to_string(ELL()));
    }
}

void Laser::buildNonzeroComponents()
{
    for (size_t idx{0}; idx < 3; ++idx)
    {
        components[idx] = (polarization[idx] != 0.0 || ell*ellipticity[idx] != 0.0) ? 1 : 0;
    }
}

double Laser::sin2_envelope(double t)
{
    double argument = (W() * t) / (2 * N());
    double value = A0() * std::sin(argument) * std::sin(argument);
    return value;
}

double Laser::A(double t, int idx)
{
    double prefactor = sin2_envelope(t) / std::sqrt(1 + ELL()*ELL());
    double term1 = Polarization()[idx] * std::sin(W()*t + CEP() - N()*M_PI);
    double term2 = ELL() * Ellipticity()[idx] * std::cos(W()*t + CEP() - N()*M_PI);

    return prefactor*(term1 + term2);
}

int Laser::Nt() const 
{
    return static_cast<int>(std::round(TMAX() / TimeSpacing())) + 1;
}

double Laser::Time(int i) const
{
    return i * TimeSpacing();
}

void Laser::dumpTo(const std::string& directory,int rank)
{
    if (rank == 0)
    {
        std::string filename = directory + "/laser.txt";

        std::ofstream outFile(filename);
    
        if (!outFile)
        {
            std::cerr << "Error opening file: " << filename << '\n';
        }
    
        for (int idx{0}; idx < Nt(); ++idx) 
        {
            outFile << Time(idx) << " " << A(Time(idx),0) << " " << A(Time(idx),1) << " " << A(Time(idx),2) << "\n";
        }
    
        outFile.close();
    }
}