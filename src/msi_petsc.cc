/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "msi_petsc.h"
#include <petsc.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <PCU.h>
#include <iostream>
namespace msi
{
  // todo : create variant using nnz structure
  msi::Mat * createPetscMatrix(int g, int l)
  {
    ::Mat * m = new ::Mat;
    MatCreateAIJ(PETSC_COMM_WORLD,
                 l,l,g,g,
                 sqrt(g), // dnz approximation
                 PETSC_NULL,
                 sqrt(g), // onz approximation
                 PETSC_NULL,
                 m);
    MatSetOption(*m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
    return reinterpret_cast<msi::Mat*>(m);
  }
  msi::Vec * createPetscVector(int g, int l)
  {
    ::Vec * v = new ::Vec;
    VecCreateMPI(PETSC_COMM_WORLD,l,g,v);
    VecSetOption(*v,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    return reinterpret_cast<msi::Vec*>(v);
  }
  ::Mat * getPetscMat(msi::Mat * m)
  {
    return reinterpret_cast<::Mat*>(m);
  }
  ::Vec * getPetscVec(msi::Vec * v)
  {
    return reinterpret_cast<::Vec*>(v);
  }
  LasOps * getPetscOps()
  {
    static PetscOps * ops = NULL;
    if(ops == NULL)
      ops = new PetscOps;
    return ops;
  }
  void PetscOps::zero(msi::Mat * m)
  {
    MatZeroEntries(*getPetscMat(m));
  }
  void PetscOps::zero(msi::Vec * v)
  {
    VecZeroEntries(*getPetscVec(v));
  }
  void PetscOps::assemble(msi::Vec * v, int cnt, int * rws, double * vls)
  {
    VecSetValues(*getPetscVec(v),cnt,rws,vls,ADD_VALUES);
  }
  void PetscOps::assemble(msi::Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls)
  {
    MatSetValues(*getPetscMat(m),cntr,rws,cntc,cls,vls,ADD_VALUES);
  }
  void PetscOps::set(msi::Vec * v, int cnt, int * rws, double * vls)
  {
    VecSetValues(*getPetscVec(v),cnt,rws,vls,INSERT_VALUES);
  }
  void PetscOps::set(msi::Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls)
  {
    MatSetValues(*getPetscMat(m),cntr,rws,cntc,cls,vls,INSERT_VALUES);
  }
  void PetscOps::multiply(msi::Mat * m, msi::Vec * x, msi::Vec * y)
  {
    MatMult(*getPetscMat(m),*getPetscVec(x),*getPetscVec(y));
  }
  double PetscOps::norm(msi::Vec * v)
  {
    double n = 0.0;
    VecNorm(*getPetscVec(v),NORM_2,&n);
    return n;
  }
  double PetscOps::dot(msi::Vec * v0, msi::Vec * v1)
  {
    double d = 0.0;
    VecDot(*getPetscVec(v0),*getPetscVec(v1),&d);
    return d;
  }
  void PetscOps::axpy(double a, Vec * x, Vec * y)
  {
    VecAXPY(*getPetscVec(y),a,*getPetscVec(x));
  }
  void PetscOps::get(msi::Vec * v, double *& vls)
  {
    VecGetArray(*getPetscVec(v),&vls);
  }
  void PetscOps::restore(msi::Vec * v, double *& vls)
  {
    VecRestoreArray(*getPetscVec(v),&vls);
  }
  class PetscLUSolve : public LasSolve
  {
  public:
    virtual void solve(msi::Mat * k, msi::Vec * u, msi::Vec * f)
    {
      ::Mat * pk = getPetscMat(k);
      ::Vec * pu = getPetscVec(u);
      ::Vec * pf = getPetscVec(f);
      ::KSP s;
      KSPCreate(PETSC_COMM_WORLD,&s);
      VecAssemblyBegin(*pf);
      VecAssemblyEnd(*pf);
      MatAssemblyBegin(*pk,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(*pk,MAT_FINAL_ASSEMBLY);
      KSPSetOperators(s,*pk,*pk);
      KSPSetFromOptions(s);
      KSPSolve(s,*pf,*pu);
      KSPDestroy(&s);
    }
    virtual int getIter() { return 0; }
  };
  LasSolve * createPetscLUSolve()
  {
    return new PetscLUSolve;
  }
  // not currently working
  class PetscQNSolve : public LasSolve
  {
  private:
    //void * args;
    int it;
  public:
    PetscQNSolve(void * a) : it(0)
      {(void) a;}
    virtual void solve(msi::Mat * k, msi::Vec * u, msi::Vec * f)
    {
      ::Mat * pk = getPetscMat(k);
      ::Vec * pu = getPetscVec(u);
      ::Vec * pf = getPetscVec(f);
      ::SNES snes;
      SNESCreate(PETSC_COMM_WORLD,&snes);
      // set snes options
      SNESSetType(snes,SNESQN);
      SNESSetJacobian(snes,*pk,*pk,NULL,NULL);
      //SNESSetFunction(snes,*pf,&PetscIterate,args);
      SNESSolve(snes,*pf,*pu);
      it = -1;
      SNESGetIterationNumber(snes,&it);
      std::cout << "BFGS QN Solve complete in " << it << " iterations." << std::endl;
      SNESDestroy(&snes);
    }
    virtual int getIter() { return it; }
  };
  LasSolve * createPetscQNSolve(void * a)
  {
    return new PetscQNSolve(a);
  }
  class PetscFieldVec : public FieldVec
  {
  protected:
    ::Vec * v;
    PetscScalar * arr;
  public:
    PetscFieldVec(apf::Numbering * n_, apf::Field * f_)
      : FieldVec(n_,f_)
      , v(NULL)
      , arr(NULL)
    {
      int l = apf::countNodes(n);
      int g = PCU_Exscan_Int(l);
      VecCreateMPI(PETSC_COMM_WORLD,l,g,v);
      VecSetOption(*v,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
      VecGetArray(*v,&arr);
      // todo: use array for field storage
      if(!apf::isFrozen(f))
        apf::freeze(f);
    }
    virtual Vec * checkout()
    {
      xd = false;
      VecRestoreArray(*v,&arr);
      return reinterpret_cast<msi::Vec*>(v);
    }
    virtual void restore(Vec *& vi)
    {
      xd = false;
      VecGetArray(*v,&arr);
      vi = NULL;
    }
  };
  FieldVec * createPetscFieldVec(apf::Numbering * nm, apf::Field * f)
  {
    return new PetscFieldVec(nm,f);
  }
}
