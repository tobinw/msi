/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifndef MSI_SOLVER_H
#define MSI_SOLVER_H
#include "msi_types.h"
namespace msi
{
  class Sim;
  Sim * createApfSim();
  void destroyApfSim(Sim * s);
  class Sim
  {
  protected:
    msi_msh * msh;
  public:
    Sim(msi_msh * m) : msh(m) {}
    virtual int registerField(msi_fld_tp tp, msi_fld * fld) = 0;
    virtual msi_fld * accessField(msi_fld_tp tp, int id) = 0;
    virtual int numberFields(msi_fld_tp tp) = 0;
    virtual msi_num * accessNumbering(msi_fld_tp tp, int id) = 0;
  };
}
#endif
