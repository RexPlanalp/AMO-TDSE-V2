#include "Angular.h"
#include <fstream>
#include <iostream>
#include "Laser.h"

void Angular::buildMaps(const Laser& laser,const std::array<int,3>& initial_state)
{   
    lm_to_block.emplace();
    block_to_lm.emplace();

    int l_i = initial_state[1];
    int m_i = initial_state[2];

    if (laser.Components()[2] && !(laser.Components()[0] || laser.Components()[1]))
    {
        buildZ(m_i);
    }
    else if ((laser.Components()[0] || laser.Components()[1]) && !laser.Components()[2])
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
