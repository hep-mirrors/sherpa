#ifndef Tensor_Current_H
#define Tensor_Current_H

#include "Vector_Current.H"
#include "C_Tensor.H"

#include <vector>

namespace EXTRAXS {

  typedef std::vector<ATOOLS::CAsT4D> CAsT_Vector;

  class Tensor_Current: public Current {
  private:
  
    CAsT_Vector m_j;

    void AddJ(const ATOOLS::CAsT4D &j);
    void Evaluate(Vertex *const v);
    void EvaluateQGC(Vector_Current *a,Vector_Current *b);

  public:

    Tensor_Current(CDBG_Amplitude *const ampl,const ATOOLS::Flavour &fl);

    ~Tensor_Current();
    
    Vertex *GetVertex(Current *const ja,Current *const jb);

    void Reset();
    void SetCurrent();
    void Evaluate();

    Complex operator*(const Current &c) const;

    void Print() const;

    inline const CAsT_Vector &J() const { return m_j; }

  };// end of class Tensor_Current

}// end of namespace EXTRAXS

#endif
