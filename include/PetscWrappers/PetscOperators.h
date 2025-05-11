#pragma once
#include "PetscWrappers/PetscVec.h"
#include "PetscWrappers/PetscMat.h"

inline Vector operator*(const Matrix& matrix, const Vector& vector)
{
    Vector temp{};
    MatCreateVecs(matrix.get(), &temp.get(),nullptr);
    MatMult(matrix.get(), vector.get(), temp.get());
    return temp;
}


inline PetscScalar innerProduct(const Vector& u, const Vector& v)
{
  PetscScalar temp{};
  VecDot(v.get(),u.get(),&temp);
  return temp;
}



inline PetscScalar norm(const Vector& vector, const Matrix& matrix)
{
  auto Sv = matrix * vector;
  PetscScalar normVal = innerProduct(Sv,vector);
  return std::sqrt(normVal);
}

inline void normalize(Vector& vector, const Matrix& matrix)
{
  auto normVal = norm(vector,matrix);
  vector *= (1.0 / normVal);
}
