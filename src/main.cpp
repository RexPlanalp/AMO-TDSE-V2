
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


#include "Simulation.h"


int main(int argc, char* argv[])
{   
    SlepcInitialize(&argc, &argv, nullptr, nullptr);
    
    MPI_Comm communicator = PETSC_COMM_WORLD;
    PetscMPIInt rank, size;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    std::string inputPath = argv[1];
    Input input{inputPath};

    TISE tise{input};
    TDSE tdse{input};
    Box box{input};
    Laser laser{input};
    Angular angular{input,laser,tdse};
    Atom atom{input};
    Basis basis{input,box};
    Observables observables{input};

    
    SimulationContext ctx = {tise,tdse,box,laser,angular,atom,basis,observables};

    Simulation simulation{size,rank,communicator,input,ctx};

    if (rank == 0)
    {   
        std::cout << tise << '\n';
        std::cout << tdse << '\n';
        std::cout << box << '\n';
        std::cout << laser << '\n';
        std::cout << angular << '\n';
        std::cout << atom << '\n';
        std::cout << basis << '\n';
        std::cout << observables << '\n';
        std::cout << simulation << '\n';

        laser.dumpTo(simulation.getMISCOutput());
        angular.dumpTo(simulation.getMISCOutput());

    }


    simulation.solveTISE();
    simulation.solveTDSE();
    simulation.computeBlockDistribution();
    simulation.computeBoundDistribution();
    simulation.computePhotoelectronSpectrum();




    
    
   


    // tdse.solve(rank,tise,basis,angular,atom,laser);
    

    //Observables observables{input};
    // observables.computeDistribution(rank,basis,tdse,tise,angular);
    // observables.computePhotoelectronSpectrum(rank,tise,tdse,angular,basis,box,atom);
    // observables.computeBoundDistribution(rank,basis,angular,tise,tdse);



    

    SlepcFinalize();
    return 0;
}