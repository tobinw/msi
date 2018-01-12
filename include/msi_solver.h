/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifndef MSI_SOLVER_H
#define MSI_SOLVER_H
#include "msi.h"
namespace msi
{
  class Sim;
  Sim * createApfSim();
  void destroyApfSim(Sim * s);
  class Sim
  {
  public:
    virtual int registerField(msi_field_type tp, msi_fld * fld) = 0;
    virtual msi_fld * accessField(msi_field_type tp, int id) = 0;
  };
}
#endif
