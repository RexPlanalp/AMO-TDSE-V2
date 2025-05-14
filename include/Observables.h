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

class Observables
{   
    public:
        Observables(const Input& input)
        : projOutBound{input.getJSON().at("Block").at("projOutBound")}
        {}
        
        void projectOutBoundStates(Vector& finalState,const Matrix& S,const TISE& tise, const Angular& angular,const Basis& Basis);

        void computeDistribution(int rank,const Basis& Basis, const TDSE& tdse,const TISE& tise, const Angular& angular);

        
       
        bool getProjOut() const {return projOutBound;}


    private:
        bool projOutBound;
};