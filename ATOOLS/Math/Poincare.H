#ifndef ATOOLS_Math_Poincare_H
#define ATOOLS_Math_Poincare_H

#include "ATOOLS/Math/Vector.H"

namespace ATOOLS {

  class Poincare {
  private:

    int m_type;
    Vec4D  m_l, m_t;
    double m_rsq, m_omct, m_st;

  public:

    Poincare(const Vec4D &v=Vec4D());
    Poincare(const Vec4D &v1,const Vec4D &v2,int mode=0);
 
    void Boost(Vec4D &v) const;
    void BoostBack(Vec4D &v) const;

    void Rotate(Vec4D &v) const;
    void RotateBack(Vec4D &v) const;

    void Lambda(Vec4D &v) const;
    void LambdaBack(Vec4D &v) const;

    void Invert();  

    inline Vec4D operator*(const Vec4D &vin) const
    {
      Vec4D v(vin);
      if (m_type==1) Boost(v);
      if (m_type==2) Rotate(v);
      if (m_type==3) Lambda(v);
      return v;
    }

    inline const Vec4D &PL() const { return m_l; }
    inline const Vec4D &PT() const { return m_t; }

    inline double OMCTheta() const { return m_omct; }
    inline double SinTheta() const { return m_st;   }

  };// end of class Poincare

  class Poincare_Sequence: public std::vector<Poincare> {
  public:

    Vec4D operator*(const Vec4D &p) const;
    void Invert();

  };// end of class Poincare_Sequence

}// end of namespace ATOOLS

#endif
