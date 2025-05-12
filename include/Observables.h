#pragma once

#include "Input.h"
#include "PetscWrappers/PetscVec.h"
#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscIS.h"
#include "PetscWrappers/PetscOperators.h"
#include "TISE.h"
#include "BSpline.h"
#include "TDSE.h"
#include "Angular.h"
#include "PetscWrappers/PetscIS.h"

class Block
{   
    public:
        Block(const Input& input)
        : projOutBound{input.getJSON().at("Block").at("projOutBound")}
        {}

        void computeDistribution(int rank,const BSpline& bspline, const TDSE& tdse,const TISE& tise, const Angular& angular);

        void projectOutBoundStates(Vector& finalState,const Matrix& S,const TISE& tise, const Angular& angular,const BSpline& bspline);
       
        bool getProjOut() const {return projOutBound;}


    private:
        bool projOutBound;
};