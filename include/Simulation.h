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
    Observables observables{};
};

struct CoulombWave
{
    std::vector<double> wave;
    double phase;
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
        , m_runBLOCK{input.getJSON().at("Simulation").at("runBlock")}
        , m_runPES{input.getJSON().at("Simulation").at("runPES")}
        , m_runBOUND{input.getJSON().at("Simulation").at("runBound")}
        , m_miscOutput{input.getJSON().at("Simulation").at("miscOutput")}
        , m_imageOutput{input.getJSON().at("Simulation").at("imageOutput")}
        , m_tiseOutput{input.getJSON().at("Simulation").at("tiseOutput")}
        , m_tdseOutput{input.getJSON().at("Simulation").at("tdseOutput")}
        {   
            buildDirectories();
        }
        
        const SimulationContext& getCtx() const {return m_ctx;}
        PetscMPIInt getRank() const {return m_rank;}
        const std::string& getTISEOutput() const {return m_tiseOutput;}
        const std::string& getTDSEOutput() const {return m_tdseOutput;}
        const std::string& getMISCOutput() const {return m_miscOutput;}
        const std::string& getIMAGEOutput() const {return m_imageOutput;}
     


        // Public Methods
        bool getTISEStatus() const {return m_runTISE;}
        bool getTDSEStatus() const {return m_runTDSE;}
        bool getHHGStatus() const {return m_runHHG;}
        bool getBLOCKStatus() const {return m_runBLOCK;}
        bool getPESStatus() const {return m_runPES;}
        bool getBOUNDStatus() const {return m_runBOUND;}

        void solveTISE();
        void solveTDSE();
        void computeBlockDistribution();
        void computeBoundDistribution();
        void computePhotoelectronSpectrum();


        


  
    private:
        // List Initialized
        PetscMPIInt m_size{};
        PetscMPIInt m_rank{};
        MPI_Comm m_communicator{};
        SimulationContext m_ctx{};
        bool m_runTISE{};
        bool m_runTDSE{};
        bool m_runHHG{};
        bool m_runBLOCK{};
        bool m_runPES{};
        bool m_runBOUND{};
        std::string m_miscOutput{};
        std::string m_imageOutput{};
        std::string m_tiseOutput{};
        std::string m_tdseOutput{};

        std::string m_eigenvalueGroup = "eigenvalues";
        std::string m_eigenvectorGroup = "eigenvectors";

        std::string m_TDSEoutputGroup = "";
        std::string m_TDSEoutputName = "psiFinal";

        // Private Methods
        void buildDirectories();

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
        void computeHHG(std::ofstream& hhgFile,int timeIdx,const Vector& state, const Matrix& XHHG, const Matrix& YHHG, const Matrix& ZHHG);

        Matrix kroneckerProduct(const Matrix& A, const Matrix& B, PetscInt nnz_A, PetscInt nnz_B);


        void projectOutBoundStates(Vector& state, const Matrix& matrix);

        CoulombWave computeCoulombWave(double E, int l);
        std::vector<std::complex<double>> expandState(const Vector& state);
        std::pair<std::map<lmPair,std::vector<std::complex<double>>>,std::map<std::pair<double, int>,double>> computePartialSpectra(const std::vector<std::complex<double>>& expanded_state);
        void computeAngleIntegrated(const std::map<lmPair,std::vector<std::complex<double>>>& partialSpectra);
        void computeAngleResolved(const std::map<lmPair,std::vector<std::complex<double>>>& partialSpectra,std::map<std::pair<double, int>,double> phases);
        
        double computeBoundPopulation(int n_bound, int l_bound, const Vector& state);
       
};   


inline std::ostream& operator<<(std::ostream& out, const Simulation& simulation)
{
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    out << "Simulation Configuration: " << "\n\n";
    out << std::setfill('\\') << std::setw(24) << "" << "\n\n";
    
    out << "TISE: " << simulation.getTISEStatus() << "\n\n";
    out << "TDSE: " << simulation.getTDSEStatus() << "\n\n";
    out << "HHG: " << simulation.getHHGStatus() << "\n\n";
    out << "BLOCK: " << simulation.getBLOCKStatus() << "\n\n";
    out << "PES: " << simulation.getPESStatus() << "\n\n";
    out << "BOUND: " << simulation.getBOUNDStatus() << "\n\n";
}
