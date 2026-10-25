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
    void Calculate(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2);
    inline double GetWeight() { return m_weight; }
    inline void SetAlphaQED(double alp) { m_alpha = alp; }
    void Subtract();
    /*!
      The W-pair Coulomb factor beyond one loop, 1 + delta_C - delta_C^(1),
      for a provider whose one-loop virtual already contains delta_C^(1)
      (flat four-fermion scheme). s: the W-pair invariant mass squared (the
      ISR-reduced s'), s1, s2: the W virtualities. order 1: O(alpha^2), eq. (7)
      of hep-ph/9507422 (Fadin, Khoze, Martin, Stirling); order 2: all orders,
      the non-relativistic factor |f|^2 of eq. (2). delta_C^(1) is eq. (9).
    */
    double HigherOrderWeight(double s, double s1, double s2, int order) const;
    //! eq. (9): the one-loop Coulomb correction delta_C^(1)
    double FirstOrder(double s, double s1, double s2) const;
  };
}

#endif
