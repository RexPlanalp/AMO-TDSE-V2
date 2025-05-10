
#include "slepc.h"

#include "Input.h"
#include "Box.h"
#include "Angular.h"
#include "Laser.h"


int main(int argc, char **argv)
{   
    
    SlepcInitialize(&argc, &argv, nullptr, nullptr);
    MPI_Comm comm = PETSC_COMM_WORLD;
    PetscMPIInt rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string inputPath = argv[1];

    // Bootstraping
    std::array<int,3> initial_state{1,0,0};
    createDirectory("misc",rank);
    createDirectory("images",rank);

    Input input{inputPath};
    input.validate();

    Box box{input};
    box.printConfiguration(rank);

    Laser laser{input};
    laser.printConfiguration(rank);
    laser.dumpTo("misc",rank);

    Angular angular{input};
    angular.buildMaps(laser,initial_state);
    angular.printConfiguration(rank);
    angular.dumpTo("misc",rank);

    SlepcFinalize();

    

    



    



    
   



  

    return 0;
}