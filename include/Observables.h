// #pragma once

// #include "Input.h"
// #include "PetscWrappers/PetscVec.h"
// #include "PetscWrappers/PetscHDF5.h"
// #include "TISE.h"
// #include "BSpline.h"

// class Block
// {   
//     public:
//         Block(const Input& input)
//         : projOutBound{input.getJSON().at("Block").at("projOutBound")}
//         {}

//         void computeDistribution(const BSpline& bspline const TDSE& tdse)
//         {
//             auto S = bspline.PopulateMatrix(PETSC_COMM_WORLD,&BSpline::overlapIntegrand, false);

//             PetscHDF5 viewer{PETSC_COMM_WORLD,}
//         }


//     private:
//         bool projOutBound;
// };