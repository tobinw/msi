/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "msi.h"
#include <pumi.h>
#include <PCU.h>
#include <parma.h>
#include <petsc.h>
#include <cassert>
#include <algorithm>
#include <functional>
#include <iostream>
bool close(const MSI_SCALAR & v1, const MSI_SCALAR & v2, const MSI_SCALAR abs_eps = 1e-8, const MSI_SCALAR rel_eps = 1e-12)
{
  if(typeid(MSI_SCALAR) == typeid(double))
  {
    double dif = abs(v1 - v2);
    bool isClose = ((dif < abs_eps) ? true : (dif <= ((v1 > v2) ? v1 : v2 * rel_eps)));
    return isClose;
  }
  else
    return false;
}
int getConfig(int argc, char * argv[], std::string & mdl_fl, std::string & msh_fl)
{
  int result = 0;
  if(argc < 3)
  {
    if(!PCU_Comm_Self())
      std::cerr << "Usage: " << argv[0] << " model(.dmg) distributed-mesh(.smb)" << std::endl;
    result = 1;
  }
  else
  {
    mdl_fl = argv[1];
    msh_fl = argv[2];
  }
  return result;
}
int main(int argc, char * argv[])
{
  MPI_Init(&argc,&argv);
  pumi_start();
  PetscInitialize(&argc,&argv,NULL,NULL);
  // suppress output to cout unless we're rank 0
  int rnk = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
  if(rnk != 0)
    std::cout.setstate(std::ios_base::failbit);
  std::string mdl_fl;
  std::string msh_fl;
  if(getConfig(argc,argv,mdl_fl,msh_fl))
  {
    // load geom and mesh however we can
    pGeom g = pumi_geom_load(mdl_fl.c_str());
    pMesh msh = pumi_mesh_load(g,msh_fl.c_str(),pumi_size());
    printStats(msh);
    int msh_dm = msi_mesh_dim(msh);
    std::cout << "* Start MSI ... " << std::endl;
    msi_start(msh);
    int dofs_per_nd = msh_dm == 3 ? 12 : 6;
    std::cout << "* creating fields - " << dofs_per_nd << " DOFs per node" << std::endl;
    msi_fld * bfld = msi_create_field(msh, "b", dofs_per_nd, 1);
    msi_fld * cfld = msi_create_field(msh, "c", dofs_per_nd, 1);
    msi_ent * lmt = msi_get_local_ent(msh,msh_dm,0);
    int nen = msi_nodes_on_ent(bfld,lmt);
    int nedofs = nen * dofs_per_nd;
    // number and create las objects

    msi_fld * xfld = msi_create_field(msh, "x", dofs_per_nd, 1);
    msi_num * num = msi_number_field(xfld,MSI_NODE_NUMBERING);
    msi_vec * xvec = msi_create_vector(num);
    msi_vec_array_field_storage(xvec,num,xfld);
    // operate on the msi_fld
    msi_vec_array_field_activate_vec(xfld);
    // operate on the msi_vec
    msi_vec_array_field_activate_field(xfld);
    // operate on the msi_fld

    
    msi_vec * bvec = msi_create_vector(num);
    msi_vec * cvec = msi_create_vector(num);
    msi_vec_as_field_storage(bvec,bfld);
    msi_vec_as_field_storage(cvec,cfld);
    msi_mat * mlt = msi_create_matrix(num);
    msi_mat * slv = msi_create_matrix(num);
    std::cout << "* set b field ..." << std::endl;
    // fill b field (and vector)
    int nv = msi_count_local_ents(msh,0);
    std::vector<MSI_SCALAR> dofs(dofs_per_nd);
    for(int idx = 0; idx < nv; ++idx)
    {
      msi_ent * vtx = msi_get_local_ent(msh,0,idx);
      int nvn = msi_nodes_on_ent(bfld,vtx);
      for(int nd = 0; nd < nvn; ++nd)
      {
        double xyz[3];
        msi_get_node_coords(bfld,vtx,nd,&xyz[0]);
        std::fill(dofs.begin(),dofs.end(),0.0);
        for(int ii = 0; ii < dofs_per_nd; ++ii)
          dofs[ii] = xyz[ii%3];
        msi_set_node_vals(bfld,vtx,nd,&dofs[0]);
      }
    }
    // fill matrices
    MSI_SCALAR ii_val = 2.0;
    MSI_SCALAR ij_val = 1.0;
    std::vector<MSI_SCALAR> blk(nedofs * nedofs,0.0);
    for(int ii = 0; ii < nedofs; ++ii)
      for(int jj = 0; jj < nedofs; ++jj)
        blk[ii*nedofs + jj] = (ii == jj ? ii_val : ij_val);
    double t1 = MPI_Wtime();
    int num_lmts = msi_count_local_ents(msh,msh_dm);
    for(int idx = 0; idx < num_lmts; ++idx)
    {
      //msi_ent * lmt = msi_get_local_ent(msh,msh_dm,idx);
      std::vector<MSI_SCALAR> tmp_blk = blk;
      for(int ii = 0; ii < 1; ++ii)
        for(int jj = 0; jj < 1; ++jj)
        {
          if(ii != jj)
            std::transform(tmp_blk.begin(), tmp_blk.end(), tmp_blk.begin(), std::bind1st(std::multiplies<MSI_SCALAR>(),0.5));
          msi_add_matrix_block(slv,ii,jj,&tmp_blk[0]);
          msi_add_matrix_block(mlt,ii,jj,&tmp_blk[0]);
        }
    }
    double t2 = MPI_Wtime();
    std::cout << "* assemble matrix ..." << std::endl;
    msi_finalize_matrix(mlt);
    msi_finalize_matrix(slv);
    double t3 = MPI_Wtime();
    std::cout << "* multiply Ab=c ..." << std::endl;
    msi_matrix_multiply(mlt, bvec, cvec);
    double t4 = MPI_Wtime();
    // field operations tests
    // copy c -> x (zero x, then x = 1.0 * c + x)
    msi_field_axpb(0.0, xfld, 0.0);
    msi_field_axpy(1.0, cfld, xfld);
    // x *= 2
    msi_field_axpb(2.0,xfld,0.0);
    // x += c
    msi_field_axpy(1.0,cfld,xfld);
    std::vector<MSI_SCALAR> cdofs(dofs_per_nd,0.0);
    std::vector<MSI_SCALAR> xdofs(dofs_per_nd,0.0);
    std::vector<MSI_SCALAR> bdofs(dofs_per_nd,0.0);
    for(int idx = 0; idx < nv; ++idx)
    {
      msi_ent * vtx = msi_get_local_ent(msh,0,idx);
      int nvn = msi_nodes_on_ent(bfld,vtx);
      for(int nd = 0; nd < nvn; ++nd)
      {
        msi_get_node_vals(xfld,vtx,nd,&xdofs[0]);
        msi_get_node_vals(cfld,vtx,nd,&cdofs[0]);
        for(int ii = 0; ii < dofs_per_nd; ++ii)
          assert(close(xdofs[ii],cdofs[ii]));
      }
    }
    // copy c field to x field
    // copy c -> x (zero x, then x = 1.0 * c + x)
    msi_field_axpb(0.0, xfld, 0.0);
    msi_field_axpy(1.0, cfld, xfld);
    double t5 = MPI_Wtime();
    std::cout << "* solve Ax=c ..." << std::endl;
    // solve Ax=c
    msi_las_solve(slv,xvec,cvec);
    double t6 = MPI_Wtime();
    // verify x=b
    std::cout << "* verify x==b..." << std::endl;
    for(int idx = 0; idx < nv; ++idx)
    {
      msi_ent * vtx = msi_get_local_ent(msh,0,idx);
      int nvn = msi_nodes_on_ent(bfld,vtx);
      for(int nd = 0; nd < nvn; ++nd)
      {
        msi_get_node_vals(xfld,vtx,nd,&xdofs[0]);
        msi_get_node_vals(bfld,vtx,nd,&bdofs[0]);
        for(int ii = 0; ii < dofs_per_nd; ++ii)
          assert(close(xdofs[ii],bdofs[ii]));
      }
    }
    std::cout << "* timings: " << std::endl
              << "  add matrix vals : " << t2-t1 << std::endl
              << "  assemble matrix : " << t3-t2 << std::endl
              << "  matrix multiply : " << t4-t3 << std::endl
              << "  solve mat sys   : " << t6-t5 << std::endl;
    msi_destroy_matrix(mlt);
    msi_destroy_matrix(slv);
    pumi_mesh_verify(msh, false);
    // should destroy the fields before the vectors used to store the fields
    msi_destroy_field(xfld);
    msi_destroy_field(bfld);
    msi_destroy_field(cfld);
    msi_destroy_vector(xvec);
    msi_destroy_vector(bvec);
    msi_destroy_vector(cvec);
    msi_finalize(msh);
    pumi_mesh_delete(msh);
  }
  pumi_finalize();
  PetscFinalize();
  MPI_Finalize();
  return 0;
}
