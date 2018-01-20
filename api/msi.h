/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
/*
  Notes: Intended for use with mdf meshes for access efficiency. Other apf and pumi
  mesh formats do not allow for direct entity access.
 */
#ifndef MSI_HEADER_H
#define MSI_HEADER_H
#include "msi_types.h"
// enums used in the API
enum msi_fld_tp { MSI_SUPPORT_FIELD = 0, MSI_SOLUTION_FIELD = 1, MSI_FIELD_TYPES = 2 };
enum msi_num_tp { MSI_DOF_NUMBERING = 0, MSI_NODE_NUMBERING = 1, MSI_NUM_TYPES = 2 };
void msi_start(msi_msh * msh);
void msi_finalize(msi_msh * msh);
// mesh api
int msi_mesh_dim(msi_msh * msh);
int msi_count_local_ents(msi_msh * msh, int dm);
msi_ent * msi_get_local_ent(msi_msh * msh, int dm, int idx);
// field api
msi_fld * msi_create_field(msi_msh * msh, const char * nm, int cmps, int ord);
void msi_destroy_field(msi_fld * fld);
msi_num * msi_number_field(msi_fld * fld, msi_num_tp tp);
int msi_count_components(msi_fld * fld);
// x = ax .+ b (element-wise)
void msi_field_axpb(MSI_SCALAR a, msi_fld * x, MSI_SCALAR b);
// y = ax + y (y != y)
void msi_field_axpy(MSI_SCALAR a, msi_fld * x, msi_fld * y);
// node api
int msi_node_id(msi_num * num, msi_ent * ent, int nd);
int msi_nodes_on_ent(msi_fld * fld, msi_ent * ent);
void msi_get_node_coords(msi_fld * fld, msi_ent * ent, int nd, double * xyz);
void msi_node_dof_range(msi_fld * fld, msi_ent * ent, int nd, int * frst, int * lst_p1);
void msi_set_node_vals(msi_fld * fld, msi_ent * ent, int nd, MSI_SCALAR * dofs);
void msi_get_node_vals(msi_fld * fld, msi_ent * ent, int nd, MSI_SCALAR * dofs);
void msi_local_node_range(msi_fld * num, int * frst, int * lst_p1);
void msi_local_node_count(msi_fld * num, int * cnt);
void msi_global_node_range(msi_fld * num, int * frst, int * lst_p1);
void msi_global_node_count(msi_fld * num, int * cnt);
// dof api
void msi_local_dof_range(msi_num * num, int * frst, int * lst_p1);
void msi_local_dof_count(msi_num * num, int * cnt);
void msi_global_dof_range(msi_num * num, int * frst, int * lst_p1);
void msi_global_dof_count(msi_num * num, int * cnt);
// las api
msi_mat * msi_create_matrix(msi_num * num);
msi_vec * msi_create_vector(msi_num * num);
void msi_destroy_matrix(msi_mat * mat);
void msi_destroy_vector(msi_vec * vec);
void msi_clear_matrix(msi_mat * mat);
void msi_clear_vector(msi_vec * vec);
void msi_set_matrix_val(msi_mat * mat, int rw, int cl, MSI_SCALAR val);
void msi_add_matrix_val(msi_mat * mat, int rw, int cl, MSI_SCALAR val);
void msi_set_vector_val(msi_vec * vec, int rw, MSI_SCALAR val);
void msi_add_vector_val(msi_vec * vec, int rw, MSI_SCALAR val);
MSI_SCALAR msi_vector_norm(msi_vec * vec);
void msi_matrix_multiply(msi_mat * a, msi_vec * x, msi_vec * y);
void msi_axpy(MSI_SCALAR * a, msi_vec * x, msi_vec * y); // y = ax+y
int msi_las_solve(msi_mat * mat, msi_vec * x, msi_vec * y);
void msi_add_matrix_block(msi_mat * mat, int brw, int bcl, MSI_SCALAR * vals);
void msi_add_matrix_blocks(msi_mat * mat, int cnt_br, int * brws, int cnt_bc, int * bcls, MSI_SCALAR * vals);
void msi_finalize_matrix(msi_mat * mat);
void msi_finalize_vector(msi_vec * vec);
// las backend storage api
void msi_vec_as_field_storage(msi_vec * vec, msi_fld * fld);
void msi_vec_array_field_storage(msi_vec * vec, msi_num * num, msi_fld * fld);
void msi_vec_array_field_activate_vec(msi_fld * fld);
void msi_vec_array_field_activate_field(msi_fld * fld);
#endif

