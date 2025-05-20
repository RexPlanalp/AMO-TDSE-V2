#include "Simulation.h"

#include "PetscWrappers/PetscHDF5.h"
#include "PetscWrappers/PetscOperators.h"

void Simulation::solveTISE()
{   
    // If we are not running the TISE, exit
    if (!getCtx().tise.getStatus())
    {
        return;
    }

    // Start measuring at the very beginning
    auto total_start = MPI_Wtime();

    // Create HDF5 File to store output
    PetscHDF5 viewer{PETSC_COMM_WORLD, getTISEOutput(), FILE_MODE_WRITE};
    
    // Create Kinetic Matrix
    auto K = Matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::K,K,false);

    auto Invr2 = Matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::Invr2,Invr2,false);
    
    auto Pot = Matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(getCtx().atom.getType(),Pot,false);

    auto S = Matrix{PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,getCtx().basis.getNbasis(),getCtx().basis.getNbasis(),2*getCtx().basis.getDegree() + 1};
    populateRadialMatrix(RadialMatrixType::S,S,false);

    K.AXPY(1.0, Pot, SAME_NONZERO_PATTERN);

    EPSSolver epssolver{PETSC_COMM_WORLD,getCtx().tise.getMaxIter(),getCtx().tise.getTol()};
    epssolver.setOperators(K,S);

    auto eigenvector = Vector{};
    S.setupVector(eigenvector);
    
    for (int l = 0; (l < getCtx().angular.getLmax()) && (l < getCtx().tise.getNmax()); ++l)
    {
        if (l > 0)
        {
            K.AXPY(l, Invr2,SAME_NONZERO_PATTERN);
        }
        
        PetscPrintf(PETSC_COMM_WORLD,"Solving for l = %d \n",l); 

        int reqPairs = getCtx().tise.getNmax() - l;
        epssolver.setDimensions(reqPairs);
        epssolver.solve();

        for (int convPair = 0; convPair < epssolver.getNconv(); ++convPair)
        {
            auto eigenvalue = epssolver.getEigenvalue(convPair);
            epssolver.getEigenvector(convPair,eigenvector);

            normalize(eigenvector,S);

            if (PetscRealPart(eigenvalue) > 0)
            {
                continue;
            }

            std::string eigenvectorName = std::string("psi_") + std::to_string(convPair + l + 1) + std::string("_") + std::to_string(l);
            std::string eigenvalueName = std::string("E_") + std::to_string(convPair + l + 1) + std::string("_") + std::to_string(l);

            viewer.saveValue(m_eigenvalueGroup, eigenvalueName, eigenvalue);
            viewer.saveVector(m_eigenvectorGroup,eigenvectorName, eigenvector);

            auto normVal = norm(eigenvector,S);
            PetscPrintf(PETSC_COMM_WORLD,"Eigenvector %d -> Norm(%.4f , %.4f) -> Eigenvalue(%.4f , %.4f)  \n",convPair+1,normVal.real(),normVal.imag(),eigenvalue.real(),eigenvalue.imag()); 
        }
    }
    auto total_end = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD,"Solved eigenvalue problem in:  %f seconds \n\n", total_end -  total_start);
}

void Simulation::populateAngularMatrix(AngularMatrixType Type, Matrix& matrix)
{ 
    for (int blockRow = 0; blockRow < getCtx().angular.getNblocks(); ++blockRow)
    {
        int l{};
        int m{};
        std::tie(l,m) = getCtx().angular.getBlockMap().at(blockRow);

        for (int blockCol = 0; blockCol < getCtx().angular.getNblocks(); ++blockCol)
        {
            int lprime{};
            int mprime{};
            std::tie(lprime,mprime) = getCtx().angular.getBlockMap().at(blockCol);

            switch(Type)
            {
                case AngularMatrixType::Z_INT_1:
                    if ((l == lprime + 1) && (m == mprime))
                    {
                        matrix.setValue(blockRow,blockCol,-PETSC_i * AngularElements::g(l,m));
                    }
                    else if ((l == lprime - 1)&&(m == mprime))
                    {   
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::f(l,m));
                    }
                break;
                case AngularMatrixType::Z_INT_2:
                    if ((l == lprime + 1) && (m == mprime))
                    {   
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::g(l,m) * (-l));
                    }
                    else if ((l == lprime - 1)&&(m == mprime))
                    {   
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::f(l,m) * (l+1));
                    }
                break;
                case AngularMatrixType::XY_INT_1:
                    if ((l == lprime + 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::a(l,m)/2);
                    }
                    else if ((l == lprime - 1)&&(m == mprime+1))
                    {   
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::b(l,m)/2);
                    }
                break;
                case AngularMatrixType::XY_INT_2:
                    if ((l == lprime + 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::c(l,m)/2);
                    }
                    else if ((l == lprime - 1)&&(m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::d(l,m)/2);
                    }
                break;
                case AngularMatrixType::XY_INT_3:
                    if ((l == lprime + 1) && (m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::atilde(l,m)/2);
                    }
                    else if ((l == lprime - 1)&&(m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::btilde(l,m)/2);
                    }
                break;
                case AngularMatrixType::XY_INT_4:
                    if ((l == lprime + 1) && (m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol, -PETSC_i * AngularElements::ctilde(l,m)/2);
                    }
                    else if ((l == lprime - 1)&&(m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol, PETSC_i * AngularElements::dtilde(l,m)/2);
                    }
                break;
                case AngularMatrixType::X_HHG:
                    if ((l == lprime + 1) && (m == mprime - 1))
                    {   
                        matrix.setValue(blockRow,blockCol,AngularElements::charlie(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime - 1))
                    {   
                        matrix.setValue(blockRow,blockCol,AngularElements::delta(l,m));
                    }
                    else if ((l == lprime + 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol,-AngularElements::alpha(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol,-AngularElements::beta(l,m));
                    }
                break;
                case AngularMatrixType::Y_HHG:
                    if ((l == lprime + 1) && (m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::charlie(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime - 1))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::delta(l,m));
                    }
                    else if ((l == lprime + 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::alpha(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime + 1))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::beta(l,m));
                    }
                break;
                case AngularMatrixType::Z_HHG:
                    if ((l == lprime + 1) && (m == mprime))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::echo(l,m));
                    }
                    else if ((l == lprime - 1) && (m == mprime))
                    {
                        matrix.setValue(blockRow,blockCol,AngularElements::foxtrot(l,m));
                    }
                break;
            }
        }
    matrix.assemble();
    }
}

void Simulation::populateRadialMatrix(RadialMatrixType Type,Matrix& matrix,bool use_ecs) 
{   
    MatrixIntegrand integrand{};

    switch(Type)
    {
        case RadialMatrixType::S:
            {
                integrand = &RadialElements::overlapIntegrand;
            }
        break;
        case RadialMatrixType::Invr2:
            {
                integrand = &RadialElements::invr2Integrand;
            }
        break;
        case RadialMatrixType::Invr:
            {
                integrand = &RadialElements::invrIntegrand;
            }
        break;
        case RadialMatrixType::Der:
            {
                integrand = &RadialElements::derIntegrand;
            }
        break;
        case RadialMatrixType::K:
            {
                integrand = &RadialElements::kineticIntegrand;
            }
        break;
        case RadialMatrixType::H:
            {
                integrand = &RadialElements::HIntegrand;
            }
        break;

    }


    for (int i = matrix.getStart(); i < matrix.getEnd(); i++) 
    {
        int col_start = std::max(0, i - getCtx().basis.getOrder() + 1);
        int col_end = std::min(getCtx().basis.getNbasis(), i + getCtx().basis.getOrder());

        for (int j = col_start; j < col_end; j++) 
        {
            std::complex<double> result = getCtx().basis.integrateMatrixElement(i, j,integrand,use_ecs);
            matrix.setValue(i,j,result);
        }
    }
    matrix.assemble();
}