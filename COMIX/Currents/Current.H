#ifndef COMIX_Currents_Current_H
#define COMIX_Currents_Current_H

#include "COMIX/Amplitude/Vertex_Base.H"
#include "METOOLS/Main/C_Scalar.H"
#include "METOOLS/Main/C_Spinor.H"
#include "METOOLS/Main/C_Vector.H"
#include "METOOLS/Main/C_Tensor.H"

namespace COMIX {

  template <typename CType> 
  class Current: public Current_Base {
  protected:

    std::vector<CType> m_j;

    friend class Current_Base;

  public:

    inline Current(const Current_Key &key):
      Current_Base(key) {}

    // inline functions
    void ResetJ();
    
    void Evaluate();

    void Print() const;

    void ResetZero();

  };// end of class Current

}// end of namespace COMIX

#endif
