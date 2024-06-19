#ifndef YFS_Main_Coulomb_H
#define YFS_Main_Coulomb_H

#include "YFS/Main/YFS_Base.H"

namespace YFS{

  class Coulomb : public YFS_Base
  {
  public:
    Coulomb();
    ~Coulomb();

    ATOOLS::Vec4D p1,p2; // final state momentum eg ttbar or W+W- 
    double m_weight, m_p, m_s1, m_s2, m_MW, m_GW, m_cms, m_threshold;
    int m_subtract;
    void Calculate(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2);
    inline double GetWeight() { return m_weight; }
    inline void SetAlphaQED(double alp) { m_alpha = alp; }
    void Subtract();
  };
}

#endif
