#include "Angular.h"
#include <fstream>
#include <iostream>
void Angular::buildMaps(const std::array<int,3>& components,const std::array<int,3>& initial_state)
{   
    lm_to_block.emplace();
    block_to_lm.emplace();

    int l_i = initial_state[1];
    int m_i = initial_state[2];

    if (components[2] && !(components[0] || components[1]))
    {
        buildZ(m_i);
    }
    else if ((components[0] || components[1]) && !components[2])
    {
        buildXY(l_i,m_i);
    }
    else
    {
        buildXYZ();
    }

    for (const auto& pair : lm_to_block.value())
    {
        block_to_lm.value()[pair.second] = pair.first;
    }

    n_lm = lm_to_block.value().size();
}

void Angular::buildZ(int m_i)
{
    for (int l = 0; l <= LMax(); ++l)
    {
        lm_to_block.value()[lm_pair(l,m_i)] = l;
    }
}

void Angular::buildXY(int l_i, int m_i)
{
    if ((l_i + m_i) % 2 == 0)
    {
        buildEven();
    }
    else
    {
        buildOdd();
    }
}

void Angular::buildXYZ()
{
    int block_idx{};

    for (int l = 0; l <= LMax(); ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            lm_to_block.value()[lm_pair(l,m)] = block_idx;
            block_idx++;
        }
    }
}

void Angular::buildEven()
{
    int block_idx{};

    for (int l = 0; l <= LMax(); ++l)
    {
        int temp_idx{};

        for (int m = -l; m <= l; ++m)
        {
            if (temp_idx % 2 == 0)
            {
                lm_to_block.value()[lm_pair(l,m)] = block_idx;
                block_idx++;
            }
            temp_idx++;
        }
    }
}

void Angular::buildOdd()
{
    int block_idx{};

    for (int l = 1; l <= LMax(); ++l)
    {
        int temp_idx{};

        for (int m = -l + 1; m <= l - 1; ++m)
        {
            if (temp_idx % 2 == 0)
            {
                lm_to_block.value()[lm_pair(l,m)] = block_idx;
                block_idx++;
            }
            temp_idx++;
        }
    }
}

void Angular::dumpTo(const std::string& directory, int rank)
{   
    if (rank == 0)
    {
        std::string filename = directory + "/lm_map.txt";

        std::ofstream outFile(filename);

        if (!outFile)
        {
            std::cerr << "Error opening file: " << filename << '\n';
        }

        for (const auto& entry : LMMap()) 
        {
            outFile << entry.first.first << " " << entry.first.second << " " << entry.second << "\n";
        }

        outFile.close();
    }
    
}

namespace AngularCoupling
{
    double f(int l, int m)
    {   
        int numerator = (l+1)*(l+1) - m*m;
        int denominator = (2*l + 1)*(2*l+3);
        return sqrt(numerator/double(denominator));
    }

    double g(int l, int m)
    {
        int numerator = l*l - m*m;
        int denominator = (2*l-1)*(2*l+1);
        return sqrt(numerator/double(denominator));
    }

    double a(int l, int m)
    {
        int numerator = (l+m);
        int denominator = (2*l +1) * (2*l-1);
        double f1 = sqrt(numerator/double(denominator));
        double f2 = - m * std::sqrt(l+m-1) - std::sqrt((l-m)*(l*(l-1)-m*(m-1)));
        return f1*f2;
        
    }

    double atilde(int l, int m)
    {
        int numerator = (l-m);
        int denominator = (2*l+1)*(2*l-1);
        double f1 = sqrt(numerator/double(denominator));
        double f2 = - m * std::sqrt(l-m-1) + std::sqrt((l+m)*(l*(l-1)-m*(m+1)));
        return f1*f2;
    }

    double b(int l, int m)
    {
        return -atilde(l+1,m-1);
    }

    double btilde(int l, int m)
    {
        return -a(l+1,m+1);
    }

    double d(int l, int m)
    {
        double numerator = (l-m+1)*(l-m+2);
        double denominator = (2*l+1)*(2*l+3);
        return std::sqrt(numerator/double(denominator));
    }

    double dtilde(int l, int m)
    {
        return d(l,-m);
    }

    double c(int l, int m)
    {
        return dtilde(l-1,m-1);
    }

    double ctilde(int l, int m)
    {
        return d(l-1,m+1);
    }

    double alpha(int l, int m)
    {
        double numerator = (l+m-1)*(l+m);
        double denominator = 4*(2*l+1)*(2*l-1);
        return std::sqrt(numerator/denominator);
    }

    double beta(int l, int m)
    {
        double numerator = (l-m+1)*(l-m+2)*(l+1);
        double denominator = 2*(2*l+1)*(2*l+2)*(2*l+3);
        return -std::sqrt(numerator/denominator);
    }

    double charlie(int l, int m)
    {
        double numerator = (l-m-1)*(l-m);
        double denominator = 4*(2*l+1)*(2*l-1);
        return std::sqrt(numerator/denominator);
    }

    double delta(int l, int m)
    {
        double numerator = (l+m+1)*(l+m+2)*(l+1);
        double denominator = 2*(2*l+1)*(2*l+2)*(2*l+3);
        return -std::sqrt(numerator/denominator);
    }

    double echo(int l, int m)
    {
        double numerator = (l+m)*(l-m);
        double denominator = (2*l-1)*(2*l+1);
        return std::sqrt(numerator/denominator);
    }

    double foxtrot(int l, int m)
    {
        double numerator = (l+m+1)*(l-m+1);
        double denominator = (2*l+1)*(2*l+3);
        return std::sqrt(numerator/denominator);
    }
}


