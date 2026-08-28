#ifndef YFS_Main_Coulomb_H
#define YFS_Main_Coulomb_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace YFS{

  class Coulomb
  {
  public:
    explicit Coulomb(bool enabled);
    ~Coulomb();

    ATOOLS::Vec4D p1,p2; // final state momentum eg ttbar or W+W- 
    double m_weight, m_p, m_s1, m_s2, m_MW, m_GW, m_cms, m_threshold;
    double m_alpha = 0.;   //!< set by SetAlphaQED, as it already was
    double m_s;            //!< sqr(Ecms), the same value YFS_Base holds
    bool   m_coulomb;
    int m_subtract;
    void Calculate(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2);
    inline double GetWeight() { return m_weight; }
    inline void SetAlphaQED(double alp) { m_alpha = alp; }
    void Subtract();
  };
}

#endif
