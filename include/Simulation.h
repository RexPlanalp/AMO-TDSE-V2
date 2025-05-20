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
        , m_runTISE{input.getJSON().at("Simulation").at("runTISE")}
        , m_runTDSE{input.getJSON().at("Simulation").at("runTDSE")}
        , m_runHHG{input.getJSON().at("Simulation").at("runHHG")}
        , m_miscOutput{input.getJSON().at("Simulation").at("miscOutput")}
        , m_imageOutput{input.getJSON().at("Simulation").at("imageOutput")}
        , m_tiseOutput{input.getJSON().at("Simulation").at("tiseOutput")}
        , m_tdseOutput{input.getJSON().at("Simulation").at("tdseOutput")}
        {
            if (getTISEStatus())
            {
                createDirectory(getTISEOutput(),getRank());
            }
        }
        
        const SimulationContext& getCtx() const {return m_ctx;}
        PetscMPIInt getRank() const {return m_rank;}
        const std::string& getTISEOutput() const {return m_tiseOutput;}


        // Public Methods
        bool getTISEStatus() const {return m_runTISE;}
        bool getTDSEStatus() const {return m_runTDSE;}
        bool getHHGStatus() const {return m_runHHG;}
        void solveTISE();
        void solveTDSE();


        


  
    private:
        // List Initialized
        PetscMPIInt m_size{};
        PetscMPIInt m_rank{};
        MPI_Comm m_communicator{};
        SimulationContext m_ctx{};
        bool m_runTISE{};
        bool m_runTDSE{};
        bool m_runHHG{};
        std::string m_miscOutput{};
        std::string m_imageOutput{};
        std::string m_tiseOutput{};
        std::string m_tdseOutput{};

        std::string m_eigenvalueGroup = "eigenvalues";
        std::string m_eigenvectorGroup = "eigenvectors";

        std::string outputGroup = "";
        std::string outputName = "psiFinal";

        // Private Methods
        void populateAngularMatrix(AngularMatrixType Type, Matrix& matrix);
        void populateRadialMatrix(RadialMatrixType Type,Matrix& matrix,bool use_ecs);

        Vector loadInitialState();
        std::pair<Matrix,Matrix> constructAtomicInteraction();
        Matrix constructZInteraction();
        std::pair<Matrix,Matrix> constructXYInteraction();
        Matrix constructXHHG();
        Matrix constructYHHG();
        Matrix constructZHHG();
        Matrix constructAtomicS();

        Matrix kroneckerProduct(const Matrix& A, const Matrix& B, PetscInt nnz_A, PetscInt nnz_B);
};      