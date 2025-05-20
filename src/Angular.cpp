#include "common.h"

#include "Angular.h"
#include "Laser.h"

void Angular::buildMaps(const Laser& laser,const TDSE& tdse)
{   
    if (laser.getComponents()[2] && !(laser.getComponents()[0] || laser.getComponents()[1]))
    {
        buildZ(tdse.getInitialM());
    }
    else if ((laser.getComponents()[0] || laser.getComponents()[1]) && !laser.getComponents()[2])
    {
        buildXY(tdse.getInitialL(),tdse.getInitialM());
    }
    else
    {
        buildXYZ();
    }

    for (const auto& pair : m_lmToBlock)
    {
        m_blockTolm[pair.second] = pair.first;
    }

    m_Nblocks = m_lmToBlock.size();
}

void Angular::buildZ(int m_i)
{
    for (int l = 0; l <= getLmax(); ++l)
    {
        m_lmToBlock[lmPair(l,m_i)] = l;
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
            m_lmToBlock[lmPair(l,m)] = block_idx;
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
                m_lmToBlock[lmPair(l,m)] = block_idx;
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
                m_lmToBlock[lmPair(l,m)] = block_idx;
                block_idx++;
            }
            temp_idx++;
        }
    }
}



// void Angular::dumpTo(const std::string& directory, int rank)
// {   
//     if (rank == 0)
//     {
//         std::string filename = directory + "/lm_map.txt";

//         std::ofstream outFile(filename);

//         if (!outFile)
//         {
//             std::cerr << "Error opening file: " << filename << '\n';
//         }

//         for (const auto& entry : getLMMap()) 
//         {
//             outFile << entry.first.first << " " << entry.first.second << " " << entry.second << "\n";
//         }

//         outFile.close();
//     }
// }
