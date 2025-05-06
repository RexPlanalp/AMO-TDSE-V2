#include <iostream>
#include <fstream>
#include <json.hpp>

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