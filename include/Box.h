#include <string>
#include "json.hpp"

class Box 
{
    public:

        Box() = delete;

        explicit Box(const nlohmann::json& input_file) : grid_size{input_file.at("Box").at("grid_size")}, grid_spacing{input_file.at("Box").at("grid_size")}
        {
            if (grid_size <= 0.0)
            {
                throw std::invalid_argument("Grid size must be greater than zero. You entered: " + std::to_string(getGridSize()));
            }
            if (grid_spacing <= 0.0)
            {
                throw std::invalid_argument("Grid spacing must be greater than zero. You entered: " + std::to_string(getGridSpacing()));
            }
        }

        double getGridSize() const {return grid_size;}
        double getGridSpacing() const {return grid_spacing;}


        int getNr() const 
        {
            return static_cast<int>(std::floor(getGridSize() / getGridSpacing()));
        }

        double getPosition(int i) const
        {
            return i * getGridSpacing();
        }

    private:
        double grid_size{};
        double grid_spacing{};


};