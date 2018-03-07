/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifndef MSI_TYPES_H_
#define MSI_TYPES_H_
//#include "msi_las.h"
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#ifdef PETSC_USE_COMPLEX
#include <complex>
struct { double dr, di; } c_omplex
typedef c_omplex MSI_SCALAR;
#else
typedef double MSI_SCALAR;
#endif
typedef apf::Mesh2 msi_msh;
typedef apf::MeshEntity msi_ent;
typedef apf::Field msi_fld;
typedef apf::FieldShape msi_shp;
typedef apf::Numbering msi_num;
namespace msi
{
  class Mat;
  class Vec;
}
typedef msi::Mat msi_mat;
typedef msi::Vec msi_vec;
#endif
