/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "msi.h"
#include "msi_field_op.h"
#include "msi_solver.h"
#include "msi_petsc.h"
#include <iostream>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <PCU.h>
#include <vector>
#include <assert.h>
#include <cmath> // isnan
msi::LasOps * ops = NULL;
//void set_adj_node_tag(pMesh m, pOwnership, pMeshTag num_global_adj_node_tag, pMeshTag num_own_adj_node_tag);
void msi_start(msi_msh * )
{
  #ifdef MSI_PETSC
  ops = msi::getPetscOps();
  #endif
}
void msi_finalize(msi_msh *)
{ }
void msi_set_node_vals(msi_fld * fld, msi_ent * ent, int nd, double * dofs)
{
  apf::setComponents(fld, ent, nd, dofs);
}
int msi_get_node_vals(msi_fld * fld, msi_ent * ent, int nd, double * dofs)
{
  apf::getComponents(fld, ent, nd, dofs);
  // need this in its own function
  return apf::countComponents(fld);
}
int mdi_node_id(msi_num * num, msi_ent * ent, int nd)
{
  return apf::getNumber(num,ent,nd,0);
}
void msi_node_dof_range(msi_num * num, msi_ent * ent, int nd, int * frst, int * lst_p1)
{
  int nd_dofs = apf::countComponents(apf::getField(num));
  int nd_id = msi_node_id(num, ent, nd);
  *frst = nd_id * nd_dofs;
  *lst_p1 = *frst + nd_dofs;
}
class PackedAXPY : public msi::FieldOp
{
protected:
  msi_fld * fld;
  msi_ent * ent;
  int cmps;
  std::vector<MSI_SCALAR> nd_dat;
  MSI_SCALAR a;
  MSI_SCALAR y;
public:
  PackedAXPY(apf::Field * f, MSI_SCALAR a_, MSI_SCALAR y_)
    : msi::FieldOp()
    , fld(f)
    , ent(NULL)
    , cmps(apf::countComponents(fld))
    , nd_dat(cmps)
    , a(a_)
    , y(y_)
  { }
  virtual bool inEntity(apf::MeshEntity * e)
  {
    ent = e;
    return true;
  }
  virtual void outEntity()
  { }
  virtual void atNode(int nd)
  {
    apf::getComponents(fld,ent,nd,&nd_dat[0]);
    for( auto val = nd_dat.begin(); val != nd_dat.end(); ++val )
    {
      *val *= a;
      *val += y;
    }
    apf::setComponents(fld,ent,nd,&nd_dat[0]);
  }
};
void msi_field_axpy(msi_fld * fld, MSI_SCALAR a, MSI_SCALAR y)
{
  PackedAXPY(fld,a,y).apply(fld);
}
msi_fld * msi_create_field(msi_msh * msh, const char * fn, int cmps, msi_shp * shp)
{
  (void)msh;
  (void)fn;
  (void)cmps;
  (void)shp;
  // todo : no direct access to tagdata
  //return apf::makeField(msh,fn,apf::PACKED,cmps,shp,new apf::TagDataOf<MSI_SCALAR>());
  return NULL;
}
void msi_local_dof_range(msi_num * num, int * frst, int * lst_p1)
{
  *lst_p1 = *frst = apf::countNodes(num) * apf::countComponents(apf::getField(num));
  *frst = PCU_Exscan_Int(*frst);
  *lst_p1 += *frst;
}
void msi_local_dof_count(msi_num * num, int * cnt)
{
  *cnt = apf::countComponents(apf::getField(num)) * apf::countNodes(num);
}
void msi_global_dof_range(msi_num * num, int * frst, int * lst_p1)
{
  *frst = 0;
  msi_global_dof_count(num,lst_p1);
}
void msi_global_dof_count(msi_num * num, int * cnt)
{
  msi_local_dof_count(num,cnt);
  *cnt = PCU_Add_Int(*cnt);
}
int msi_count_components(msi_fld * fld)
{
  return apf::countComponents(fld);
}
msi_mat * msi_create_matrix(msi_num * num)
{
  msi_mat * mat = NULL;
  int gbl = 0;
  int lcl = 0;
  msi_global_dof_count(num,&gbl);
  msi_local_dof_count(num,&lcl);
  #ifdef USE_PETSC
  mat = msi::createPetscMatrix(gbl,lcl);
  #endif
  return mat;
}
void msi_destroy_matrix(msi_mat * mat)
{
  msi::destroyPetscMatrix(mat);
}
void msi_destroy_vector(msi_vec * vec)
{
  msi::destroyPetscVector(vec);
}
void msi_clear_matrix(msi_mat * mat)
{
  ops->zero(mat);
}
void msi_clear_vector(msi_vec * vec)
{
  ops->zero(vec);
}
void msi_set_matrix_val(msi_mat * mat, int rw, int cl, MSI_SCALAR val)
{
  ops->set(mat,1,&rw,1,&cl,&val);
}
void msi_add_matrix_val(msi_mat * mat, int rw, int cl, MSI_SCALAR val)
{
  ops->assemble(mat,1,&rw,1,&cl,&val);
}
void msi_set_vector_val(msi_vec * vec, int rw, MSI_SCALAR val)
{
  ops->set(vec,1,&rw,&val);
}
void msi_add_vector_val(msi_vec * vec, int rw, MSI_SCALAR val)
{
  ops->assemble(vec,1,&rw,&val);
}
void msi_matrix_multiply(msi_mat * mat, msi_vec * x, msi_vec * y)
{
  ops->multiply(mat,x,y);
}
void msi_axpy(MSI_SCALAR * a, msi_vec * x, msi_vec * y)
{
  ops->axpy(*a,x,y);
}
int msi_las_solve(msi_mat * mat, msi_vec * x, msi_vec * y)
{
  msi::LasSolve * slvr = msi::createPetscLUSolve();
  slvr->solve(mat,x,y);
  return slvr->getIter();
}
void msi_add_matrix_block(msi_mat * mat, msi_num  * num, msi_ent * ent, MSI_SCALAR * vals)
{
  apf::NewArray<int> dof_ids;
  int dof_cnt = apf::getElementNumbers(num, ent, dof_ids);
  ops->assemble(mat,dof_cnt,&dof_ids[0],dof_cnt,&dof_ids[0],vals);
}
/*
struct entMsg
{
  int pid;
  apf::MeshEntity* ent;
  entMsg( int pid_p=0, apf::MeshEntity* ent_p=NULL)
  {
    pid=pid_p;
    ent=ent_p;
  }
};
struct classcomp
{
  bool operator() (const entMsg& lhs, const entMsg& rhs) const
  {
    if(lhs.ent==rhs.ent) return lhs.pid<rhs.pid;
    else return lhs.ent<rhs.ent;
  }
};
void set_adj_node_tag(pMesh m, pOwnership o, pMeshTag num_global_adj_node_tag, pMeshTag num_own_adj_node_tag)
{
  int value;
  int brgType = m->getDimension()-1;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  PCU_Comm_Begin();
  while ((e = m->iterate(it)))
  {
    int num_adj_node=0;
    apf::Adjacent elements;
    apf::getBridgeAdjacent(m, e, brgType, 0, elements);
    int num_adj = elements.getSize();
    for (int i=0; i<num_adj; ++i)
    {
      if (pumi_ment_isOwned(elements[i], o))
        ++num_adj_node;
    }
    m->setIntTag(e, num_own_adj_node_tag, &num_adj_node);
    if (!m->isShared(e)) continue;
    // first pass msg size to owner
    int own_partid = pumi_ment_getOwnPID(e,o);
    apf::MeshEntity* own_copy = pumi_ment_getOwnEnt(e,o);
    if (!own_copy) // own_copy does not exist so let;'s
    {
    }
    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_copy);
    PCU_Comm_Pack(own_partid, &num_adj,sizeof(int));
  }
  m->end(it);
  PCU_Comm_Send();
  std::map<apf::MeshEntity*, std::map<int, int> > count_map;
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(e);
      PCU_Comm_Unpack(&value,sizeof(int));
      count_map[e][PCU_Comm_Sender()]=value;
    }
  }
  // pass entities to owner
  std::map<apf::MeshEntity*, std::set<entMsg, classcomp> > count_map2;
  it = m->begin(0);
  PCU_Comm_Begin();
  while ((e = m->iterate(it)))
  {
    // pass entities to ownner
    std::vector<entMsg> msgs;
    apf::Adjacent elements;
    apf::getBridgeAdjacent(m, e, brgType, 0, elements);
    apf::MeshEntity* ownerEnt=pumi_ment_getOwnEnt(e,o);
    int own_partid = pumi_ment_getOwnPID(e, o);
    for(int i=0; i<elements.getSize(); ++i)
    {
      apf::MeshEntity* ownerEnt2=pumi_ment_getOwnEnt(elements[i],o);
      int owner=pumi_ment_getOwnPID(elements[i], o);
      msgs.push_back(entMsg(owner, ownerEnt2));
      if(own_partid==PCU_Comm_Self())
      {
        count_map2[e].insert(*msgs.rbegin());
      }
    }
    if(own_partid!=PCU_Comm_Self())
    {
      PCU_COMM_PACK(own_partid, ownerEnt);
      PCU_Comm_Pack(own_partid, &msgs.at(0),sizeof(entMsg)*msgs.size());
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(e);
      int sizeData = count_map[e][PCU_Comm_Sender()];
      std::vector<entMsg> data(sizeData);
      PCU_Comm_Unpack(&data.at(0),sizeof(entMsg)*sizeData);
      for (int i=0; i<data.size(); ++i)
      {
        count_map2[e].insert(data.at(i));
      }
    }
  }
  for (std::map<apf::MeshEntity*, std::set<entMsg,classcomp> >::iterator mit=count_map2.begin();
       mit!=count_map2.end(); ++mit)
  {
    e = mit->first;
    int num_global_adj =count_map2[e].size();
    m->setIntTag(mit->first, num_global_adj_node_tag, &num_global_adj);
  }
}
*/
