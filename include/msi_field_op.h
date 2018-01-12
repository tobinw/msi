/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#ifndef MSI_FIELD_OP_H_
#define MSI_FIELD_OP_H_
#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
namespace msi
{
  class FieldOp
  {
  public:
    virtual bool inEntity(apf::MeshEntity *) {return true;}
    virtual void outEntity() {}
    virtual void atNode(int) {}
    void apply(apf::Field * f)
    {
      apf::Mesh * msh = apf::getMesh(f);
      apf::FieldShape * s = apf::getShape(f);
      for(int ii = 0; ii < 4; ++ii)
      {
        if(!s->hasNodesIn(ii))
          continue;
        apf::MeshIterator * it = msh->begin(ii);
        apf::MeshEntity * e = NULL;
        while((e = msh->iterate(it)))
        {
          int n = s->countNodesOn(msh->getType(e));
          if(n==0)
            continue;
          if(!this->inEntity(e))
            continue;
          for(int jj = 0; jj < n; ++jj)
            this->atNode(jj);
          this->outEntity();
        }
        msh->end(it);
      }
    }
  };
}
#endif
