#pragma once

#include <string>
#include "mpi.h"
#include "nlohmann/json.hpp"

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

        explicit Simulation(const std::string& inputPath, int ranks, MPI_Comm comm)
        : inputParams{loadJson(inputPath)}, comm{comm}, ranks{ranks}
        , box{inputParams}
        , angular{inputParams}
        , bspline{inputParams,box}
        , atom{inputParams}
        , laser{inputParams}
        {}
    
    private:

        Box box;
        Angular angular;
        BSpline bspline;
        Laser laser;
        Atom atom;

        nlohmann::json inputParams{};
        MPI_Comm comm{};
        int ranks{};
};