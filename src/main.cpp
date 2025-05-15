
#include "slepc.h"

#include "Input.h"
#include "Box.h"
#include "Angular.h"
#include "Laser.h"
#include "Atom.h"
#include "Basis.h"
#include "TISE.h"
#include "TDSE.h"
#include "Observables.h"


#include <chrono>
#include "PetscWrappers/PetscMat.h"
#include "petsc.h"


int main(int argc, char **argv)
{   
    
    SlepcInitialize(&argc, &argv, nullptr, nullptr);
    MPI_Comm comm = PETSC_COMM_WORLD;
    PetscMPIInt rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string inputPath = argv[1];

    
    
    Input input{rank,inputPath};

    Box box{input};
    box.printConfiguration(rank);

    Laser laser{input};
    laser.printConfiguration(rank);
    laser.dumpTo("misc",rank);

    Angular angular{input};
    angular.printConfiguration(rank);
    angular.dumpTo("misc",rank);

    Atom atom{input};
    atom.printConfiguration(rank);

    Basis basis{input};
    basis.buildKnots(box);
    basis.printConfiguration(rank);
    //basis,dumpTo(box,"misc",rank);

    TISE tise{rank,input};
    tise.printConfiguration(rank);
    tise.solve(basis,atom,angular);

   

    TDSE tdse{rank,input};
    tdse.printConfiguration(rank);

    angular.buildMaps(laser,tdse.getInitialNLM());


    tdse.solve(rank,tise,basis,angular,atom,laser);
    

    Observables observables{input};
    observables.computeDistribution(rank,basis,tdse,tise,angular);
    observables.computePhotoelectronSpectrum(rank,tise,tdse,angular,basis,box,atom);



    

    SlepcFinalize();
    return 0;
}