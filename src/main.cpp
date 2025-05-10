
#include "slepc.h"

#include "Input.h"
#include "Box.h"


int main(int argc, char **argv)
{   
    
    SlepcInitialize(&argc, &argv, nullptr, nullptr);
    MPI_Comm comm = PETSC_COMM_WORLD;
    PetscMPIInt rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string inputPath = argv[1];

    Input input{inputPath};
    input.validate();

    Box box{input};
    box.printConfiguration(rank);

    



    
   



  

    return 0;
}