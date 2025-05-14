#include "common.h"

#include "Angular.h"
#include "Laser.h"

void Angular::buildMaps(const Laser& laser,const std::array<int,3>& initial_state)
{   
    lm_to_block.emplace();
    block_to_lm.emplace();

    int l_i = initial_state[1];
    int m_i = initial_state[2];

    if (laser.getComponents()[2] && !(laser.getComponents()[0] || laser.getComponents()[1]))
    {
        buildZ(m_i);
    }
    else if ((laser.getComponents()[0] || laser.getComponents()[1]) && !laser.getComponents()[2])
    {
        buildXY(l_i,m_i);
    }
    else
    {
        buildXYZ();
    }

    for (const auto& pair : lm_to_block)
    {
        block_to_lm[pair.second] = pair.first;
    }

    nlm = lm_to_block.size();
}

void Angular::buildZ(int m_i)
{
    for (int l = 0; l <= getLmax(); ++l)
    {
        lm_to_block[lm_pair(l,m_i)] = l;
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

    for (int l = 0; l <= getLmax(); ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            lm_to_block[lm_pair(l,m)] = block_idx;
            block_idx++;
        }
    }
}

void Angular::buildEven()
{
    int block_idx{};

    for (int l = 0; l <= getLmax(); ++l)
    {
        int temp_idx{};

        for (int m = -l; m <= l; ++m)
        {
            if (temp_idx % 2 == 0)
            {
                lm_to_block[lm_pair(l,m)] = block_idx;
                block_idx++;
            }
            temp_idx++;
        }
    }
}

void Angular::buildOdd()
{
    int block_idx{};

    for (int l = 1; l <= getLmax(); ++l)
    {
        int temp_idx{};

        for (int m = -l + 1; m <= l - 1; ++m)
        {
            if (temp_idx % 2 == 0)
            {
                lm_to_block[lm_pair(l,m)] = block_idx;
                block_idx++;
            }
            temp_idx++;
        }
    }
}

void Angular::printConfiguration(int rank)
{
    if (rank == 0)
    {   
        std::cout << std::setfill('\\') << std::setw(24) << "" << '\n';
        std::cout << "Angular Configuration: " << "\n\n";
        std::cout << std::setfill('\\') << std::setw(24) << "" << '\n';

        std::cout << "lmax: " << getLmax() << "\n\n";
        std::cout << "mmax: " << getMmax() << "\n\n";
        std::cout << "mmin: " << getMmin() << "\n\n";
        std::cout << "nlm: "  << getNlm() << "\n\n";
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

        for (const auto& entry : getLMMap()) 
        {
            outFile << entry.first.first << " " << entry.first.second << " " << entry.second << "\n";
        }

        outFile.close();
    }
}
