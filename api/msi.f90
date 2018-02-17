module msi
  interface
#ifdef USECOMPLEX
#define MSI_SCALAR complex(c_double_complex)
#else
#define MSI_SCALAR real(c_double)
#endif
     !
     !  Util api
     !
     subroutine msi_start() bind(c)
       use iso_c_binding
     end subroutine msi_start
     subroutine msi_finalize() bind(c)
       use iso_c_binding
     end subroutine msi_finalize
     !
     !  Mesh api
     !
     integer(c_int) function msi_mesh_dim(msh) bind(c)
       use iso_c_binding
       type(c_ptr), value :: msh
     end function msi_mesh_dim
     integer(c_int) function msi_count_local_ents(msh,dm) bind(c)
       use iso_c_binding
       type(c_ptr), value :: msh
       integer(c_int), value :: dm
     end function msi_count_local_ents
     type(c_ptr) function msi_get_local_ent(msh,dm,idx) bind(c)
       use iso_c_binding
       type(c_ptr), value :: msh
       integer(c_int), value :: dm
       integer(c_int), value :: idx
     end function msi_get_local_ent
     integer(c_int) function msi_get_ent_id(msh,ent) bind(c)
       use iso_c_binding
       type(c_ptr), value :: msh
       type(c_ptr), value :: ent
     end function msi_get_ent_id
     subroutine msi_get_adjacent(msh,ent,adj_dim,adj) bind(c)
       use iso_c_binding
       type(c_ptr), value :: msh
       type(c_ptr), value :: ent
       integer(c_int), value :: adj_dim
       type(c_ptr)        :: adj
     end subroutine msi_get_adjacent
     subroutine msi_get_geom_class(msh,ent,mdl_dim,mdl_id) bind(c)
       use iso_c_binding
       type(c_ptr), value :: msh
       type(c_ptr), value :: ent
       integer(c_int), intent(out) :: mdl_dim
       integer(c_int), intent(out) :: mdl_id
     end subroutine msi_get_geom_class
     !
     !  Field api
     !
     type(c_ptr) function msi_create_field(msh,nm,cmps,ord) bind(c)
       use iso_c_binding
       type(c_ptr), value :: msh
       char(c_char), intent(in) :: nm(:)
       integer(c_int), value :: cmps
       integer(c_int), value :: ord
     end function msi_create_field
     subroutine msi_destroy_field(fld) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
     end subroutine msi_destroy_field
     type(c_ptr) function msi_number_field(fld,tp) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
       integer(c_int), value :: tp
     end function msi_number_field
     integer(c_int) function msi_count_components(fld) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
     end function msi_count_components
     subroutine msi_field_axpb(a,x,b) bind(c)
       use iso_c_binding
       MSI_SCALAR, value :: a
       type(c_ptr), value :: x
       MSI_SCALAR, value :: b
     end subroutine msi_field_axpb
     subroutine msi_field_axpy(a,x,y) bind(c)
       use iso_c_binding
       MSI_SCALAR, value :: a
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end subroutine msi_field_axpy
     !
     ! Node api
     !
     integer(c_int) function msi_node_id(num,ent,nd) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
       type(c_ptr), value :: ent
       integer(c_int), value :: nd
     end function msi_node_id
     integer(c_int) function msi_nodes_on_ent(fld,ent) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
       type(c_ptr), value :: ent
     end function msi_nodes_on_ent
     subroutine msi_get_node_coords(fld,ent,nd,xyz) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
       type(c_ptr), value :: ent
       integer(c_int), value :: nd
       real(c_double), intent(out) :: xyz(3)
     end subroutine msi_get_node_coords
     subroutine msi_node_dof_range(fld,ent,nd,frst,lst_p1) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
       type(c_ptr), value :: ent
       integer(c_int), value :: nd
       integer(c_int), intent(out) :: frst
       integer(c_int), intent(out) :: lst_p1
     end subroutine msi_node_dof_range
     subroutine msi_set_node_vals(fld,ent,nd,dofs) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
       type(c_ptr), value :: ent
       integer(c_int), value :: nd
       MSI_SCALAR, intent(in) :: dofs(:)
     end subroutine msi_set_node_vals
     subroutine msi_get_node_vals(fld,ent,nd,dofs) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
       type(c_ptr), value :: ent
       integer(c_int), value :: nd
       MSI_SCALAR, intent(out) :: dofs(:)
     end subroutine msi_get_node_vals
     subroutine msi_local_node_range(num,frst,lst_p1) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
       integer(c_int), intent(out) :: frst
       integer(c_int), intent(out) :: lst_p1
     end subroutine msi_local_node_range
     integer(c_int) function msi_local_node_count(num) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
     end function msi_local_node_count
     subroutine msi_global_node_range(num,frst,lst_p1) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
       integer(c_int), intent(out) :: frst
       integer(c_int), intent(out) :: lst_p1
     end subroutine msi_global_node_range
     integer(c_int) function msi_global_node_count(num) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
     end function msi_global_node_count
     !
     !  Dof api
     !
     subroutine msi_local_dof_range(num,frst,lst_p1) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
       integer(c_int), intent(out) :: frst
       integer(c_int), intent(out) :: lst_p1
     end subroutine msi_local_dof_range
     integer(c_int) function msi_local_dof_count(num) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
     end function msi_local_dof_count
     subroutine msi_global_dof_range(num,frst,lst_p1) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
       integer(c_int), intent(out) :: frst
       integer(c_int), intent(out) :: lst_p1
     end subroutine msi_global_dof_range
     integer(c_int) function msi_global_dof_count(num) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
     end function msi_global_dof_count
     !
     ! Las api
     !
     type(c_ptr) function msi_create_matrix(num) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
     end function msi_create_matrix
     type(c_ptr) function msi_create_vector(num) bind(c)
       use iso_c_binding
       type(c_ptr), value :: num
     end function msi_create_vector
     subroutine msi_destroy_matrix(mat) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
     end subroutine msi_destroy_matrix
     subroutine msi_destroy_vector(vec) bind(c)
       use iso_c_binding
       type(c_ptr), value :: vec
     end subroutine msi_destroy_vector
     subroutine msi_clear_matrix(mat) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
     end subroutine msi_clear_matrix
     subroutine msi_clear_vector(vec) bind(c)
       use iso_c_binding
       type(c_ptr), value :: vec
     end subroutine msi_clear_vector
     subroutine msi_set_matrix_val(mat,rw,cl,val) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       integer(c_int), value :: rw
       integer(c_int), value :: cl
       MSI_SCALAR, value :: val
     end subroutine msi_set_matrix_val
     subroutine msi_add_matrix_val(mat,rw,cl,val) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       integer(c_int), value :: rw
       integer(c_int), value :: cl
       MSI_SCALAR, value :: val
     end subroutine msi_add_matrix_val
     subroutine msi_set_vector_val(vec,rw,val) bind(c)
       use iso_c_binding
       type(c_ptr), value :: vec
       integer(c_int), value :: rw
       MSI_SCALAR, value :: val
     end subroutine msi_set_vector_val
     subroutine msi_add_vector_val(vec,rw,val) bind(c)
       use iso_c_binding
       type(c_ptr), value :: vec
       integer(c_int), value :: rw
       MSI_SCALAR, value :: val
     end subroutine msi_add_vector_val
     MSI_SCALAR function msi_vector_norm(vec) bind(c)
       use iso_c_binding
       type(c_ptr), value :: vec
     end function msi_vector_norm
     subroutine msi_matrix_multiply(a,x,y) bind(c)
       use iso_c_binding
       type(c_ptr), value :: a
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end subroutine msi_matrix_multiply
     subroutine msi_axpy(a,x,y) bind(c)
       use iso_c_binding
       MSI_SCALAR, value :: a
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end subroutine msi_axpy
     integer(c_int) function msi_las_solve(mat,x,y) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end function msi_las_solve
     subroutine msi_add_matrix_block(mat,brw,bcl,vals) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       integer(c_int), value :: brw
       integer(c_int), value :: bcl
       MSI_SCALAR, intent(in) :: vals
     end subroutine msi_add_matrix_block
     subroutine msi_add_matrix_blocks(mat,brw_cnt,brws,bcl_cnt,bcls,vals) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       integer(c_int), value :: brw_cnt
       integer(c_int), intent(in) :: brws
       integer(c_int), value :: bcl_cnt
       integer(c_int), intent(in) :: bcls
       MSI_SCALAR, intent(in) :: vals
     end subroutine msi_add_matrix_blocks
     subroutine msi_finalize_matrix(mat) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
     end subroutine msi_finalize_matrix
     subroutine msi_finalize_vector(vec) bind(c)
       use iso_c_binding
       type(c_ptr), value :: vec
     end subroutine msi_finalize_vector
     !
     ! Las storage api
     !
     subroutine msi_vec_as_field_storage(vec,fld) bind(c)
       use iso_c_binding
       type(c_ptr), value :: vec
       type(c_ptr), value :: fld
     end subroutine msi_vec_as_field_storage
     subroutine msi_vec_array_field_storage(vec,num,fld) bind(c)
       use iso_c_binding
       type(c_ptr), value :: vec
       type(c_ptr), value :: num
       type(c_ptr), value :: fld
     end subroutine msi_vec_array_field_storage
     subroutine msi_vec_array_field_activate_vec(fld) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
     end subroutine msi_vec_array_field_activate_vec
     subroutine msi_vec_array_field_activate_field(fld) bind(c)
       use iso_c_binding
       type(c_ptr), value :: fld
     end subroutine msi_vec_array_field_activate_field
  end interface
contains
end module msi
