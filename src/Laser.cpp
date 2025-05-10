#include "Laser.h"


void Laser::buildNonzeroComponents()
{
    for (size_t idx{0}; idx < 3; ++idx)
    {
        components[idx] = (polarization[idx] != 0.0 || ell*ellipticity[idx] != 0.0) ? 1 : 0;
    }
}

double Laser::sin2_envelope(double t)
{
    double argument = (getW() * t) / (2 * getN());
    double value = getA0() * std::sin(argument) * std::sin(argument);
    return value;
}

double Laser::A(double t, int idx)
{
    double prefactor = sin2_envelope(t) / std::sqrt(1 + getEll()*getEll());
    double term1 = getPolarization()[idx] * std::sin(getW()*t + getCEP() - getN()*M_PI);
    double term2 = getEll() * getEllipticity()[idx] * std::cos(getW()*t + getCEP() - getN()*M_PI);

    return prefactor*(term1 + term2);
}

double Laser::getTime(int i) const
{
    return times[i];
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
    
        for (int idx{0}; idx < getNt(); ++idx) 
        {
            outFile << getTime(idx) << " " << A(getTime(idx),0) << " " << A(getTime(idx),1) << " " << A(getTime(idx),2) << "\n";
        }
    
        outFile.close();
    }
}

void Laser::printConfiguration(int rank)
{
    if (rank == 0)
    {
        std::cout << "Laser Configuration: " << "\n\n";
        std::cout << "N_cycles: " << getN() << "\n\n";
        std::cout << "dt: " << getTimeSpacing() << "\n\n";
        std::cout << "tmax: " << getTmax() << "\n\n";
        std::cout << "Nt: " << getNt() << "\n\n";
        std::cout << "w: " << getW() << "\n\n";
        std::cout << "I: " << getI() << "\n\n";
        std::cout << "A_0: " << getA0() << "\n\n";
        std::cout << "ell: " << getEll() << "\n\n";
        std::cout << "CEP: " << getCEP() << "\n\n";
        std::cout << "polarization: " << getPolarization()[0] << " " << getPolarization()[1] << " " << getPolarization()[2] <<  "\n\n";
        std::cout << "poynting: " << getPoynting()[0] << " " << getPoynting()[1] << " " << getPoynting()[2] <<  "\n\n";
        std::cout << "ellipticity: " << getEllipticity()[0] << " " << getEllipticity()[1] << " " << getEllipticity()[2] <<  "\n\n";
        std::cout << "components: " << getComponents()[0] << " " << getComponents()[1] << " " << getComponents()[2] <<  "\n\n";


    }
}