#pragma once

#include "petscksp.h"
#include "PetscWrappers/PetscMat.h"
#include "PetscWrappers/PetscVec.h"

class KSPSolver
{
    public:
        KSPSolver(MPI_Comm comm, PetscInt maxIter, PetscReal tolerance, PetscInt restart)
        : maxIter{maxIter}
        , tolerance{tolerance}
        {
            KSPCreate(comm,&ksp);
            KSPSetType(ksp, KSPGMRES);
            KSPSetTolerances(ksp,tolerance,PETSC_DETERMINE,PETSC_DETERMINE,maxIter);
            KSPGMRESSetRestart(ksp, restart);
            KSPSetReusePreconditioner(ksp,PETSC_TRUE);
            KSPSetFromOptions(ksp);
        }

        ~KSPSolver()
        {
            if (ksp)
            {
                KSPDestroy(&ksp);
            }
        }

        void setOperators(const Matrix& left)
        {
            KSPSetOperators(ksp,left.get(),left.get());
        }

        void setPreconditioner(PetscInt Nblocks)
        {
            PC pc;
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCBJACOBI);  
            PCBJacobiSetTotalBlocks(pc, Nblocks, nullptr);
        }

        void solve(Vector& rhs)
        {
            KSPSolve(ksp, rhs.get(),rhs.get());
        }

        const KSP& get() const {return ksp;}
        KSP& get() {return ksp;}
       

    private:
        PetscInt maxIter{};
        PetscReal tolerance{};

        PetscInt restart{};


        KSP ksp{nullptr};
};