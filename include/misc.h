#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>

inline nlohmann::json loadJson(const std::string& path)
{
    std::ifstream file(path);

    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open json file: " + path);
    }

    nlohmann::json j{};

    try 
    {
        file >> j;
    }
    catch(const std::exception& e)
    {
        throw std::runtime_error("Failed to parse json file: " + path + "\n" + e.what());
    }

    return j;
}

inline void createDirectory(int rank, const std::string& directory)
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
