#include "Laser.h"

void Laser::nonzeroComponents()
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