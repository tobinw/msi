/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "msi.h"
#include "msi_solver.h"
#include <apf.h>
#include <apfShape.h>
#include <apfMesh.h>
#include <vector>
namespace msi
{
  class apf_sim : public Sim
  {
  public:
    virtual int registerField(msi_fld_tp tp, msi_fld * fld)
    {
      flds[tp].push_back(fld);
      return static_cast<int>(flds[tp].size());
    }
    virtual msi_fld * accessField(msi_fld_tp tp, int id)
    {
      // assert id is < size
      return flds[tp].at(id);
    }
  protected:
    std::vector<apf::Field*> flds[MSI_FIELD_TYPES];
  };
}
