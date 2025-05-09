#pragma once

#include <string>
#include "mpi.h"

#include "Box.h"
#include "Angular.h"
#include "BSpline.h"
#include "Laser.h"
#include "Atom.h"
#include "misc.h"

class Simulation
{   
    public:
        Simulation() = delete;

        explicit Simulation(const std::string& inputPath, int current_rank, int ranks,  MPI_Comm comm)
        : inputParams{loadJson(inputPath)}, comm{comm}, ranks{ranks}, current_rank{current_rank}
        , box{loadJson(inputPath)}
        , angular{loadJson(inputPath)}
        , bspline{loadJson(inputPath),box}
        , atom{loadJson(inputPath)}
        , laser{loadJson(inputPath)}
        
        {

            std::array<int,3> initial_state= loadJson(inputPath).at("TDSE").at("initial_state").get<std::array<int,3>>();
            angular.buildMaps(laser,initial_state);
        }
    
       
        
    
    private:
        nlohmann::json inputParams{};
        MPI_Comm comm{};
        int ranks{};
        int current_rank{};
    
    public:
        Box box;
        Angular angular;
        BSpline bspline;
        Atom atom;
        Laser laser;
};      