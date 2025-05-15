#pragma once

#include "Input.h"
#include "PetscWrappers/PetscVec.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscIS.h"
#include "PetscWrappers/PetscOperators.h"
#include "TISE.h"
#include "Basis.h"
#include "TDSE.h"
#include "Angular.h"
#include "PetscWrappers/PetscIS.h"

struct CoulombWave
{
    std::vector<double> wave;
    double phase;
};

class Observables
{   
    public:
        Observables(const Input& input)
        : projOutBound{input.getJSON().at("Observables").at("Block").at("projOutBound")}
        , Emin{input.getJSON().at("Observables").at("PES").at("Emin")}
        , Emax{input.getJSON().at("Observables").at("PES").at("Emax")}
        , slice{input.getJSON().at("Observables").at("PES").at("slice")}
        , pesStatus{input.getJSON().at("Observables").at("PES").at("status")}
        , blockStatus{input.getJSON().at("Observables").at("Block").at("status")}
        , boundStatus{input.getJSON().at("Observables").at("Bound").at("status")}
        {
             Ne = static_cast<int>(std::round(getEmax() / getEmin())) + 1;
        }
        
        void projectOutBoundStates(Vector& finalState,const Matrix& S,const TISE& tise, const Angular& angular,const Basis& basis);

        void computeDistribution(int rank,const Basis& basis, const TDSE& tdse,const TISE& tise, const Angular& angular);

        CoulombWave computeCoulombWave(double E, int l, const Box& box, const Atom& atom);
        std::vector<std::complex<double>> expandState(const Vector& state, const Box& box, const Angular& angular, const Basis& basis);
        std::pair<std::map<lm_pair,std::vector<std::complex<double>>>,std::map<std::pair<double, int>,double>> computePartialSpectra(const std::vector<std::complex<double>>& expanded_state,const Angular& angular, const Atom& atom,const Box& box);
        void computeAngleIntegrated(const std::map<lm_pair,std::vector<std::complex<double>>>& partialSpectra,const Angular& angular);
        void computeAngleResolved(const std::map<lm_pair,std::vector<std::complex<double>>>& partialSpectra,std::map<std::pair<double, int>,double> phases);
        void computePhotoelectronSpectrum(int rank,const TISE& tise, const TDSE& tdse, const Angular& angular, const Basis& basis, const Box& box, const Atom& atom);

        double computeBoundPopulation(int n_bound, int l_bound, const Vector& state,const TISE& tise,const Basis& basis,const Angular& angular);
        void computeBoundDistribution(int rank,const Basis& basis, const Angular& angular, const TISE& tise, const TDSE& tdse);
        
       
        bool getProjOut() const {return projOutBound;}
        double getEmin() const {return Emin;}
        double getEmax() const {return Emax;}
        const std::string& getSlice() const {return slice;}
        int getNe() const {return Ne;}
        bool getPESStatus() const {return pesStatus;}
        bool getBlockStatus() const {return blockStatus;}
        bool getBoundStatus() const {return boundStatus;}

        void printConfiguration(int rank);


    private:
        bool projOutBound{};
        double Emin{};
        double Emax{};
        int Ne{};
        std::string slice{};

        bool pesStatus{};
        bool blockStatus{};
        bool boundStatus{};
};
