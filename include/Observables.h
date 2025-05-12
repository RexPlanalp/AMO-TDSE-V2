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

class Block
{   
    public:
        Block(const Input& input)
        : projOutBound{input.getJSON().at("Block").at("projOutBound")}
        {}

        void computeDistribution(int rank,const BSpline& bspline, const TDSE& tdse, const Angular& angular);
       


    private:
        bool projOutBound;
};