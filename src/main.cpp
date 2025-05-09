

#include "Simulation.h"
#include "misc.h"
#include "Box.h"
#include "Angular.h"
#include "Laser.h"
#include "BSpline.h"
#include "Atom.h"

#include <array>
#include <nlohmann/json.hpp>
#include "Simulation.h"

#include "slepc.h"


int main(int argc, char **argv)
{   
    
    SlepcInitialize(&argc, &argv, nullptr, nullptr);
    MPI_Comm comm = PETSC_COMM_WORLD;
    PetscMPIInt rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string inputPath = "input.json";

    Simulation simulation(inputPath, rank, size, comm);

   

    simulation.bspline.dumpTo(simulation.box,"misc",rank);
    simulation.angular.dumpTo("misc",rank);
    simulation.laser.dumpTo("misc",rank);

  

    return 0;
}