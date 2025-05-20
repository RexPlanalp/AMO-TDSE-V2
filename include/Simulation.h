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

        explicit Simulation(PetscMPIInt size, PetscMPIInt rank, MPI_Comm communicator, const Input& input, const SimulationContext& ctx)
        : m_size{size}
        , m_rank{rank}
        , m_communicator{communicator}
        , m_ctx{ctx}
        , m_miscOutput{input.getJSON().at("Simulation").at("miscOutput")}
        , m_imageOutput{input.getJSON().at("Simulation").at("imageOutput")}
        , m_tiseOutput{input.getJSON().at("Simulation").at("tiseOutput")}
        , m_tdseOutput{input.getJSON().at("Simulation").at("tdseOutput")}
        {
            if (getCtx().tise.getStatus())
            {
                createDirectory(getTISEOutput(),getRank());
            }
        }
        
        const SimulationContext& getCtx() const {return m_ctx;}
        PetscMPIInt getRank() const {return m_rank;}
        const std::string& getTISEOutput() const {return m_tiseOutput;}


        // Public Methods
        void solveTISE();

  
    private:
        // List Initialized
        PetscMPIInt m_size{};
        PetscMPIInt m_rank{};
        MPI_Comm m_communicator{};
        SimulationContext m_ctx{};
        std::string m_miscOutput{};
        std::string m_imageOutput{};
        std::string m_tiseOutput{};
        std::string m_tdseOutput{};

        std::string m_eigenvalueGroup = "eigenvalues";
        std::string m_eigenvectorGroup = "eigenvectors";

        // Private Methods
        void populateAngularMatrix(AngularMatrixType Type, Matrix& matrix);
        void populateRadialMatrix(RadialMatrixType Type,Matrix& matrix,bool use_ecs);
};      