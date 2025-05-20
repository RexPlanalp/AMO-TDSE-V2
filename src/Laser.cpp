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

void Laser::dumpTo(const std::string& directory)
{

    std::string filename = directory + "/laser.txt";

    std::ofstream outFile(filename);

    if (!outFile)
    {
        std::cerr << "Error opening file: " << filename << '\n';
    }

    for (int idx{0}; idx < getNt(); ++idx) 
    {
        outFile << (*this)[idx] << " " << A((*this)[idx],0) << " " << A((*this)[idx],1) << " " << A((*this)[idx],2) << "\n";
    }

    outFile.close();

}
