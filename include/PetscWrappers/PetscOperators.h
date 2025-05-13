#pragma once
#include "PetscWrappers/PetscVec.h"
#include "PetscWrappers/PetscMat.h"



inline void matMult(const Matrix& m, const Vector& input, Vector& output)
{
  MatMult(m.get(),input.get(),output.get());
}

inline PetscScalar innerProduct(const Vector& u, const Vector& v)
{
  PetscScalar temp{};
  VecDot(v.get(),u.get(),&temp);
  return temp;
}

inline PetscScalar innerProduct(const Vector& u, const Matrix& m, const Vector& v)
{

  auto temp = Vector{};
  m.setupVector(temp);

  m.matMult(v,temp);


  return innerProduct(u, temp);
}



inline PetscScalar norm(const Vector& vector, const Matrix& matrix)
{
  PetscScalar normVal = innerProduct(vector,matrix,vector);
  return std::sqrt(normVal);
}

inline void normalize(Vector& vector, const Matrix& matrix)
{
  auto normVal = norm(vector,matrix);
  vector *= (1.0 / normVal);
}
