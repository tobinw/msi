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
#include <cassert>
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
  inline ::Mat * getPetscMat(msi::Mat * m)
  {
    return reinterpret_cast<::Mat*>(m);
  }
  inline ::Vec * getPetscVec(msi::Vec * v)
  {
    return reinterpret_cast<::Vec*>(v);
  }
  void destroyPetscMatrix(msi::Mat * m)
  {
    MatDestroy(getPetscMat(m));
  }
  void destroyPetscVector(msi::Vec * v)
  {
    VecDestroy(getPetscVec(v));
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
    assert(x != y);
    ::Mat * pm = getPetscMat(m);
    ::Vec * vx = getPetscVec(x);
    ::Vec * vy = getPetscVec(y);
    PetscBool ass;
    MatAssembled(*pm,&ass);
    if(!ass)
    {
      MatAssemblyBegin(*pm,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(*pm,MAT_FINAL_ASSEMBLY);
    }
    VecAssemblyBegin(*vx);
    VecAssemblyEnd(*vx);
    VecAssemblyBegin(*vy);
    VecAssemblyEnd(*vy);
    MatMult(*pm,*vx,*vy);
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
  void PetscOps::finalize(msi::Mat * m)
  {
    ::Mat * pm = getPetscMat(m);
    PetscBool ass;
    MatAssembled(*pm,&ass);
    if(!ass)
    {
      MatAssemblyBegin(*pm,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(*pm,MAT_FINAL_ASSEMBLY);
    }
  }
  void PetscOps::finalize(msi::Vec * v)
  {
    ::Vec * pv = getPetscVec(v);
    VecAssemblyBegin(*pv);
    VecAssemblyEnd(*pv);
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
      PetscBool ass;
      MatAssembled(*pk,&ass);
      if(!ass)
      {
        MatAssemblyBegin(*pk,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*pk,MAT_FINAL_ASSEMBLY);
      }
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
  template <class T>
  class PetscVecDataOf : public apf::FieldDataOf<T>
  {
  protected:
    msi_vec * vec;
    msi_num * num;
  public:
    PetscVecDataOf(msi_vec * v)
      : apf::FieldDataOf<T>()
      , vec(v)
      , num(NULL)
    { }
    virtual void init(apf::FieldBase * f)
    {
      this->field = f;
      apf::FieldShape * s = f->getShape();
      const char * name = s->getName();
      apf::Numbering * n = f->getMesh()->findNumbering(name);
      if (!n)
        n = numberOverlapNodes(f->getMesh(),name,s);
      num = n;
    }
    virtual ~PetscVecDataOf() { }
    virtual bool hasEntity(apf::MeshEntity*)
    {
      return true;
    }
    virtual void removeEntity(apf::MeshEntity*)
    { }
    virtual void get(apf::MeshEntity * e, T * data)
    {
      apf::NewArray<int> dofs;
      int dof_cnt = apf::getElementNumbers(num,e,dofs);
      VecGetValues(*getPetscVec(vec),dof_cnt,&dofs[0],&data[0]);
    }
    virtual void set(apf::MeshEntity * e, T const * data)
    {
      apf::NewArray<int> dofs;
      int dof_cnt = apf::getElementNumbers(num,e,dofs);
      VecSetValues(*getPetscVec(vec),dof_cnt,&dofs[0],&data[0],INSERT_VALUES);
    }
    virtual bool isFrozen()
    {
      return true;
    }
  };
  template <class T>
  class PetscVecArrayDataOf : public apf::FieldDataOf<T>
  {
  protected:
    msi_vec * vec;
    T * arr;
    msi_num * num;
    int frst;
    bool valid;
  public:
    PetscVecArrayDataOf(msi_vec * v, msi_num * n)
      : apf::FieldDataOf<T>()
      , vec(v)
      , arr(NULL)
      , num(n)
      , frst(0)
      , valid(true)
    { }
    virtual void init(apf::FieldBase * f)
    {
      this->field = f;
      frst = PCU_Exscan_Int(apf::countNodes(num) * apf::countComponents(apf::getField(num)));
      // assert that the vector has enough size to store the data
    }
    virtual ~PetscVecArrayDataOf() { }
    virtual bool hasEntity(apf::MeshEntity * ent)
    {
      return apf::getMesh(num)->isOwned(ent);
    }
    virtual void removeEntity(apf::MeshEntity*)
    { }
    virtual void get(apf::MeshEntity * e, T * data)
    {
      if(!valid)
        std::cerr << "Error! Cannot access field data while underlying vector is active! Call .activateField()" << std::endl;
      apf::NewArray<int> dofs;
      int dof_cnt = apf::getElementNumbers(num,e,dofs);
      // need to offset the dof vals by the first local dof value...
      for(int ii = 0; ii < dof_cnt; ++ii)
        data[ii] = arr[dofs[ii]-frst];
    }
    virtual void set(apf::MeshEntity * e, T const * data)
    {
      if(!valid)
        std::cerr << "Error! Cannot access field data while underlying vector is active! Call .activateField()" << std::endl;
      apf::NewArray<int> dofs;
      int dof_cnt = apf::getElementNumbers(num,e,dofs);
      // could replace with memcpy if all accesses are contiguous
      for(int ii = 0; ii < dof_cnt; ++ii)
        arr[dofs[ii]-frst] = data[ii];
    }
    virtual bool isFrozen()
    {
      return true;
    }
    void activateVec()
    {
      VecRestoreArray(*getPetscVec(vec),&arr);
      valid = false;
    }
    void activateField()
    {
      VecGetArray(*getPetscVec(vec),&arr);
      valid = true;
    }
  };
  template <class T>
  void vectorFieldData(msi_fld * fld, msi_vec * vec)
  {
    PetscVecDataOf<T>* newData = new PetscVecDataOf<T>(vec);
    newData->init(fld);
    apf::FieldDataOf<T> * oldData = static_cast<apf::FieldDataOf<T>*>(fld->getData());
    apf::copyFieldData<T>(oldData,newData);
    fld->changeData(newData);
  }
  template <class T>
  void vectorArrayFieldData(msi_fld * fld, msi_num * num, Vec * vec)
  {
    PetscVecArrayDataOf<T>* newData = new PetscVecArrayDataOf<T>(vec,num);
    newData->init(fld);
    apf::FieldDataOf<T> * oldData = static_cast<apf::FieldDataOf<T>*>(fld->getData());
    apf::copyFieldData<T>(oldData,newData);
    fld->changeData(newData);
  }
  template <class T>
  void vectorArrayFieldActivateField(msi_fld * fld)
  {
    PetscVecArrayDataOf<T>* data = static_cast<PetscVecArrayDataOf<T>*>(fld->getData());
    data->activateField();
  }
  template <class T>
  void vectorArrayFieldActivateVec(msi_fld * fld)
  {
    PetscVecArrayDataOf<T>* data = static_cast<PetscVecArrayDataOf<T>*>(fld->getData());
    data->activateVec();
  }
  //instantiate the template functions
  template void vectorFieldData<MSI_SCALAR>(msi_fld * fld, msi_vec * vec);
  template void vectorArrayFieldData<MSI_SCALAR>(msi_fld * fld, msi_num * num, msi_vec * vec);
  template void vectorArrayFieldActivateField<MSI_SCALAR>(msi_fld * fld);
  template void vectorArrayFieldActivateVec<MSI_SCALAR>(msi_fld * fld);
}
