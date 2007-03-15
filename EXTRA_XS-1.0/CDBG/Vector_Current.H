#ifndef Vector_Current_H
#define Vector_Current_H

#include "Current.H"
#include "C_Vector.H"

#include <vector>

namespace EXTRAXS {

  typedef std::vector<ATOOLS::CVec4D> CVec_Vector;

  class Tensor_Current;

  class Vector_Current: public Current {
  private:
  
    CVec_Vector m_j;

    void AddJ(const ATOOLS::CVec4D &j);
    void Evaluate(Vertex *const v);
    void EvaluateTGC(Vector_Current *a,Vector_Current *b);
    void EvaluateQGC(Vector_Current *a,Tensor_Current *b,
		     const int mode);

  public:

    Vector_Current(CDBG_Amplitude *const ampl,const ATOOLS::Flavour &fl);

    ~Vector_Current();
    
    Vertex *GetVertex(Current *const ja,Current *const jb);

    void Reset();
    void SetCurrent();
    void Evaluate();

    Complex operator*(const Current &c) const;

    void Print() const;

    inline const CVec_Vector &J() const { return m_j; }

  };// end of class Vector_Current

}// end of namespace EXTRAXS

#endif
