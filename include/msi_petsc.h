/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifndef MSI_PETSC_H
#define MSI_PETSC_H
#include "msi.h"
#include "msi_solver.h"
#include "msi_field_storage.h"
#include "msi_las.h"
#include <apf.h>
#include <apfNumbering.h>
#include <vector>
#include <map>
namespace msi
{
  Mat * createPetscMatrix(int gbl, int lcl);
  void destroyPetscMatrix(Mat * m);
  Vec * createPetscVector(int gbl, int lcl);
  void destroyPetscVector(Vec * v);
  FieldVec * createPetscFieldVec(apf::Numbering * nm, apf::Field * f);
  LasOps * getPetscOps();
  LasSolve * createPetscLUSolve();
  LasSolve * createPetscQNSolve(void * a);
  class PetscOps : public LasOps
  {
  public:
    virtual void zero(Mat * m);
    virtual void zero(Vec * v);
    virtual void assemble(Vec * v, int cnt, int * rws, double * vls);
    virtual void assemble(Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls);
    virtual void set(Vec * v, int cnt, int * rws, double * vls);
    virtual void set(Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls);
    virtual void multiply(Mat * m, Vec * x, Vec * y);
    virtual double norm(Vec * v);
    virtual double dot(Vec * v0, Vec * v1);
    virtual void axpy(double a, Vec * x, Vec * y);
    virtual void get(Vec * v, double *& vls);
    virtual void restore(Vec * v, double *& vls);
  };
}
#endif
