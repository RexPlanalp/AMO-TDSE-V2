#pragma once

#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <complex>


namespace Constants
{
    constexpr double I_AU = 3.51E16;
}

inline nlohmann::json loadJson(const std::string& path)
{
    auto abs = std::filesystem::absolute(path);
    std::ifstream file(abs);

    if (!file.is_open()) 
    {
      throw std::runtime_error("…");
    }

    nlohmann::json j;
    file >> j;
    return j;
}

inline void createDirectory(const std::string& directory,int rank)
{
    try 
    {
        std::filesystem::create_directories(directory);
    } catch (const std::filesystem::filesystem_error& e) 
    {
        throw std::runtime_error(
            "Rank " + std::to_string(rank)
          + ": Failed to create directory '" + directory
          + "': " + e.what()
        );
    }
}

template<typename T>
inline void normalize(std::array<T,3>& array)
{
    T norm = std::sqrt(array[0] * array[0] + array[1] * array[1] + array[2] * array[2]);

    if (norm <= T{})
    {
        throw std::domain_error("Cannot divide by zero.");
    }

    for (size_t idx{0}; idx < array.size(); ++idx)
    {
        array[idx] /= norm;
    }
}

template<typename T>
inline std::array<T,3> crossProduct(const std::array<T,3>& array1, const std::array<T,3>& array2)
{   
    std::array<T,3> result{};
    result[0] = array1[1]*array2[2] - array1[2]*array2[1];
    result[1] = array1[2]*array2[0] - array1[0]*array2[2];
    result[2] = array1[0]*array2[1] - array1[1]*array2[0];

    return result;
}

template<typename T> 
inline double realDotProduct(const T& container1, const T& container2)
{
    if (container1.size() != container2.size())
    {
        throw std::invalid_argument("Containers must be the same length.");
    }
    

    double innerProduct{};

    for (size_t idx{0}; idx < container1.size(); ++idx)
    {
        innerProduct += container1[idx] * container2[idx];
    }

    return innerProduct;
}


inline std::complex<double> SimpsonsMethod(const std::vector<std::complex<double>>& vec,double dr)
{
    int n = vec.size() - 1;
    std::complex<double> I {};


    if (n % 2 != 0)
    {
        I += (3 * dr / 8) * (vec[n-2]+static_cast<double>(3)*vec[n-1]+static_cast<double>(3)*vec[n]);
        n -= 2;
    }

    double p = dr / 3; 
    I += (vec[0]+vec[n]) * p;
    for (int vec_idx = 1; vec_idx < n; vec_idx += 2)
    {
        I += static_cast<double>(4) * vec[vec_idx] * p;
    }
    for (int vec_idx = 2; vec_idx < n-1; vec_idx += 2)
    {
        I += static_cast<double>(2) * vec[vec_idx] * p;
    }
    
    return I;

}