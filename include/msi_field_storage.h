/******************************************************************************
  (c) 2017 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
*******************************************************************************/
#include "msi_las.h"
#include <apf.h>
#include <apfNumbering.h>
namespace msi
{
// use the underlying storage from a vector to store a field, this class receives a field, creates a parallel vector of sufficient size to store the field, copies the field data out, and replaces the underlying field data array with the underlying array from the vector
  class FieldVec
  {
  public:
    FieldVec(apf::Numbering * n_, apf::Field * f_)
      : n(n_)
      , f(f_)
      , xd(true)
    { }
    virtual Vec * checkout() = 0;
    virtual void restore(Vec *&) = 0;
    virtual bool checkedOut()
    {
      return xd;
    }
  protected:
    apf::Numbering * n;
    apf::Field * f;
    bool xd;
  };
}
