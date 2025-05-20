#include "Laser.h"


void Laser::buildParameters()
{
    m_I /= Constants::I_AU;
    m_A_0 = std::sqrt(m_I) / getW();

    m_Tmax = getN() * 2 * M_PI / getW();
    m_cep *= M_PI;
}

void Laser::buildNt()
{
    m_Nt = static_cast<int>(std::round(getTmax() / getTimeSpacing())) + 1;
}

void Laser::buildTimes()
{
    m_times.resize(getNt());
    for (int idx = 0; idx < getNt(); ++idx)
    {
        m_times[idx] = idx * getTimeSpacing();
    }
}

void Laser::buildVectors()
{
    m_ellipticity = crossProduct(getPolarization(),getPoynting());
    normalize(m_polarization);
    normalize(m_poynting);
    normalize(m_ellipticity);

    for (size_t idx{0}; idx < 3; ++idx)
    {
        m_components[idx] = (getPolarization()[idx] != 0.0 || getEll()*getEllipticity()[idx] != 0.0) ? 1 : 0;
    }
}


double Laser::sin2_envelope(double t) const
{
    double argument = (getW() * t) / (2 * getN());
    double value = getA0() * std::sin(argument) * std::sin(argument);
    return value;
}

double Laser::A(double t, int idx) const
{
    double prefactor = sin2_envelope(t) / std::sqrt(1 + getEll()*getEll());
    double term1 = getPolarization()[idx] * std::sin(getW()*t + getCEP() - getN()*M_PI);
    double term2 = getEll() * getEllipticity()[idx] * std::cos(getW()*t + getCEP() - getN()*M_PI);

    return prefactor*(term1 + term2);
}

double Laser::operator[](int i) const
{
    return m_times[i];
}

// void Laser::dumpTo(const std::string& directory,int rank)
// {
//     if (rank == 0)
//     {
//         std::string filename = directory + "/laser.txt";

//         std::ofstream outFile(filename);
    
//         if (!outFile)
//         {
//             std::cerr << "Error opening file: " << filename << '\n';
//         }
    
//         for (int idx{0}; idx < getNt(); ++idx) 
//         {
//             outFile << getTime(idx) << " " << A(getTime(idx),0) << " " << A(getTime(idx),1) << " " << A(getTime(idx),2) << "\n";
//         }
    
//         outFile.close();
//     }
// }

// void Laser::printConfiguration(int rank)
// {
//     if (rank == 0)
//     {
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
//         std::cout << "Laser Configuration: " << "\n\n";
//         std::cout << std::setfill('\\') << std::setw(24) << "" << "\n\n";
        
//         std::cout << "N_cycles: " << getN() << "\n\n";
//         std::cout << "dt: " << getTimeSpacing() << "\n\n";
//         std::cout << "tmax: " << getTmax() << "\n\n";
//         std::cout << "Nt: " << getNt() << "\n\n";
//         std::cout << "w: " << getW() << "\n\n";
//         std::cout << "I: " << getI() << "\n\n";
//         std::cout << "A_0: " << getA0() << "\n\n";
//         std::cout << "ell: " << getEll() << "\n\n";
//         std::cout << "CEP: " << getCEP() << "\n\n";
//         std::cout << "polarization: " << getPolarization()[0] << " " << getPolarization()[1] << " " << getPolarization()[2] <<  "\n\n";
//         std::cout << "poynting: " << getPoynting()[0] << " " << getPoynting()[1] << " " << getPoynting()[2] <<  "\n\n";
//         std::cout << "ellipticity: " << getEllipticity()[0] << " " << getEllipticity()[1] << " " << getEllipticity()[2] <<  "\n\n";
//         std::cout << "components: " << getComponents()[0] << " " << getComponents()[1] << " " << getComponents()[2] <<  "\n\n";


//     }
// }