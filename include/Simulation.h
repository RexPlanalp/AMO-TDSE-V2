#pragma once

#include "nlohmann/json.hpp"
#include <string>
#include "mpi.h"

#include "Box.h"
#include "Angular.h"
#include "BSpline.h"
#include "Laser.h"
#include "Atom.h"
#include "misc.h"
#include "TISE.h"

class Simulation
{   
    public:
        Simulation() = delete;

        explicit Simulation(const nlohmann::json& inputPar, int rank, int size,  MPI_Comm comm)
        : comm{comm}, size{size}, rank{rank}
        , box{inputPar}
        , angular{inputPar}
        , bspline{inputPar,box}
        , atom{inputPar}
        , laser{inputPar}
        , tise{inputPar}
        
        {

            createDirectory(rank, "misc");
            createDirectory(rank, "images");

            std::array<int,3> initial_state= inputPar.at("TDSE").at("initial_state").get<std::array<int,3>>();
            angular.buildMaps(laser,initial_state);

            bspline.dumpTo(box,"misc",rank);
            angular.dumpTo("misc",rank);
            laser.dumpTo("misc",rank);
        }

        void solveTISE()
        {
            tise.solve(bspline,atom, angular);
        }
    
       
        
    
    private:
        // Member List Initialized
        MPI_Comm comm{};
        int size{};
        int rank{};

        // Default InitializedBox box;
        Box box;
        Angular angular;
        BSpline bspline;
        Atom atom;
        Laser laser;
        TISE tise;
};      