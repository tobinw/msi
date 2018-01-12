/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifndef MSI_HEADER_H
#define MSI_HEADER_H
#include "msi_types.h"
void msi_start(msi_msh * msh);
void msi_finalize(msi_msh * msh);
// Field api
msi_fld * msi_create_field(msi_msh * msh, const char * nm, int comp);
int msi_count_components(msi_fld * fld);
void msi_field_axpy(msi_fld * fld, MSI_SCALAR a, MSI_SCALAR y);
// Node api
int msi_node_id(msi_num * num, msi_ent * ent, int nd);
void msi_node_dof_range(msi_fld * fld, msi_ent * ent, int nd, int * frst, int * lst_p1);
void msi_set_node_vals(msi_fld * fld, msi_ent * ent, int nd, double * dofs);
int msi_get_node_vals(msi_fld * fld, msi_ent * ent, int nd, double * dofs);
// dof api
void msi_local_dof_range(msi_num * num, int * frst, int * lst_p1);
void msi_local_dof_count(msi_num * num, int * cnt);
void msi_global_dof_range(msi_num * num, int * frst, int * lst_p1);
void msi_global_dof_count(msi_num * num, int * cnt);
// las api
msi_mat * msi_create_matrix(msi_fld * num);
msi_vec * msi_create_vector(msi_fld * num);
void msi_destroy_matrix(msi_mat * mat);
void msi_destory_vector(msi_vec * vec);
void msi_clear_matrix(msi_mat * mat);
void msi_clear_vector(msi_vec * vec);
void msi_set_matrix_val(msi_mat * mat, int rw, int cl, MSI_SCALAR val);
void msi_add_matrix_val(msi_mat * mat, int rw, int cl, MSI_SCALAR val);
void msi_set_vector_val(msi_vec * vec, int rw, MSI_SCALAR val);
void msi_add_vector_val(msi_vec * vec, int rw, MSI_SCALAR val);
void msi_matrix_multiply(msi_mat * a, msi_vec * x, msi_vec * y);
// y = ax+y
void msi_axpy(MSI_SCALAR * a, msi_vec * x, msi_vec * y);
int msi_las_solve(msi_mat * mat, msi_fld * x, msi_fld * y);
void msi_add_matrix_block(msi_mat * mat, msi_num * num, msi_ent * ent, int rwidx, int clidx, MSI_SCALAR * vals);
#endif

