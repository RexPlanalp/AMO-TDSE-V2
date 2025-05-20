#pragma once

#include "common.h"

#include "Input.h"


struct SimulationContext
{   
    TISE tise{};
    TDSE tdse{};
    Box box{};
    Laser laser{};
    Angular angular{};
    Atom atom{};
    Basis basis{};
};


class Simulation
{   
    public:

        explicit Simulation(PetscMPIInt size, PetscMPIInt rank, MPI_Comm communicator, SimulationContext ctx)
        : m_size{size}
        , m_rank{rank}
        , m_communicator{communicator}
        , m_ctx{ctx}
        {}
  
        
        

       
        
    
    private:
        PetscMPIInt m_size{};
        PetscMPIInt m_rank{};
        MPI_Comm m_communicator{};

        SimulationContext m_ctx{};

      
};      