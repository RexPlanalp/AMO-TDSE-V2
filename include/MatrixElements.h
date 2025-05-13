#pragma once

#include "common.h"
#include "BSplines.h"
#include "PetscWrappers/PetscMat.h"
#include "Angular.h"
#include "Basis.h"

enum class AngularMatrixType
{
    Z_INT_1,
    Z_INT_2,
    XY_INT_1,
    XY_INT_2,
    XY_INT_3,
    XY_INT_4
};

enum class RadialMatrixType
{
    S,
    Invr,
    Invr2,
    K,
    Der,
    H
};

namespace AngularElements
{
    inline double f(int l, int m)
    {   
        int numerator = (l+1)*(l+1) - m*m;
        int denominator = (2*l + 1)*(2*l+3);
        return std::sqrt(numerator/double(denominator));
    }

    inline double g(int l, int m)
    {
        int numerator = l*l - m*m;
        int denominator = (2*l-1)*(2*l+1);
        return std::sqrt(numerator/double(denominator));
    }

    inline double a(int l, int m)
    {
        int numerator = (l+m);
        int denominator = (2*l +1) * (2*l-1);
        double f1 = sqrt(numerator/double(denominator));
        double f2 = - m * std::sqrt(l+m-1) - std::sqrt((l-m)*(l*(l-1)-m*(m-1)));
        return f1*f2;
        
    }

    inline double atilde(int l, int m)
    {
        int numerator = (l-m);
        int denominator = (2*l+1)*(2*l-1);
        double f1 = sqrt(numerator/double(denominator));
        double f2 = - m * std::sqrt(l-m-1) + std::sqrt((l+m)*(l*(l-1)-m*(m+1)));
        return f1*f2;
    }

    inline double b(int l, int m)
    {
        return -atilde(l+1,m-1);
    }

    inline double btilde(int l, int m)
    {
        return -a(l+1,m+1);
    }

    inline double d(int l, int m)
    {
        double numerator = (l-m+1)*(l-m+2);
        double denominator = (2*l+1)*(2*l+3);
        return std::sqrt(numerator/double(denominator));
    }

    inline double dtilde(int l, int m)
    {
        return d(l,-m);
    }

    inline double c(int l, int m)
    {
        return dtilde(l-1,m-1);
    }

    inline double ctilde(int l, int m)
    {
        return d(l-1,m+1);
    }

    inline double alpha(int l, int m)
    {
        double numerator = (l+m-1)*(l+m);
        double denominator = 4*(2*l+1)*(2*l-1);
        return std::sqrt(numerator/denominator);
    }

    inline double beta(int l, int m)
    {
        double numerator = (l-m+1)*(l-m+2)*(l+1);
        double denominator = 2*(2*l+1)*(2*l+2)*(2*l+3);
        return -std::sqrt(numerator/denominator);
    }

    inline double charlie(int l, int m)
    {
        double numerator = (l-m-1)*(l-m);
        double denominator = 4*(2*l+1)*(2*l-1);
        return std::sqrt(numerator/denominator);
    }

    inline double delta(int l, int m)
    {
        double numerator = (l+m+1)*(l+m+2)*(l+1);
        double denominator = 2*(2*l+1)*(2*l+2)*(2*l+3);
        return -std::sqrt(numerator/denominator);
    }

    inline double echo(int l, int m)
    {
        double numerator = (l+m)*(l-m);
        double denominator = (2*l-1)*(2*l+1);
        return std::sqrt(numerator/denominator);
    }

    inline double foxtrot(int l, int m)
    {
        double numerator = (l+m+1)*(l-m+1);
        double denominator = (2*l+1)*(2*l+3);
        return std::sqrt(numerator/denominator);
    }
}

namespace RadialElements
{
    inline std::complex<double> overlapIntegrand(int i, int j, int degree, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) {return BSplines::B(i,degree, x,localKnots) * BSplines::B(j,degree, x,localKnots);}
    inline std::complex<double> kineticIntegrand(int i, int j,int degree, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) {return 0.5 * BSplines::dB(i,degree, x,localKnots) * BSplines::dB(j,degree, x,localKnots);}
    inline std::complex<double> invrIntegrand(int i, int j,int degree, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) {return BSplines::B(i,degree, x,localKnots) * BSplines::B(j,degree, x,localKnots) / (x + 1E-25);}
    inline std::complex<double> invr2Integrand(int i, int j,int degree, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) {return BSplines::B(i,degree, x,localKnots) * BSplines::B(j,degree, x,localKnots) / (x*x + 1E-25);}
    inline std::complex<double> derIntegrand(int i, int j,int degree, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) {return BSplines::B(i,degree, x,localKnots) * BSplines::dB(j,degree,x,localKnots);}
    inline std::complex<double> HIntegrand(int i, int j,int degree, std::complex<double> x, const std::vector<std::complex<double>>& localKnots) {return  BSplines::B(i,degree, x,localKnots) * BSplines::B(j,degree, x,localKnots) * Potentials::hydrogen(x) ;}
}

namespace AngularMatrix
{
    inline void populateAngularMatrix(AngularMatrixType Type, Matrix& matrix, const Angular& angular)
    { 
        for (int blockRow = 0; blockRow < angular.getNlm(); ++blockRow)
        {
            int l{};
            int m{};
            std::tie(l,m) = angular.getBlockMap().at(blockRow);

            for (int blockCol = 0; blockCol < angular.getNlm(); ++blockCol)
            {
                int lprime{};
                int mprime{};
                std::tie(lprime,mprime) = angular.getBlockMap().at(blockCol);

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
                }
            }
        matrix.assemble();
        }
    }
}

namespace RadialMatrix
{
    
    inline void populateRadialMatrix(RadialMatrixType Type,Matrix& matrix,const Basis& basis,bool use_ecs) 
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
            int col_start = std::max(0, i - basis.getOrder() + 1);
            int col_end = std::min(basis.getNbasis(), i + basis.getOrder());

            for (int j = col_start; j < col_end; j++) 
            {
                std::complex<double> result = basis.integrateMatrixElement(i, j,integrand,use_ecs);
                matrix.setValue(i,j,result);
            }
        }
        matrix.assemble();
    }
}