

#include "Simulation.h"
#include "misc.h"
#include "Box.h"
#include "Angular.h"
#include "Laser.h"
#include "BSpline.h"
#include "Atom.h"

#include <array>
#include "misc.h"
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

    std::string inputPath = argv[1];

    nlohmann::json inputPar = loadJson(inputPath);

    Simulation simulation(inputPar, rank, size, comm);
    simulation.solveTISE();
   



  

    return 0;
}