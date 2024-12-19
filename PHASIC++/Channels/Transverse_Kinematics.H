#ifndef PHASIC__Channels__Transverse_Kinematics_H
#define PHASIC__Channels__Transverse_Kinematics_H

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  struct LN_Pair {
    ATOOLS::Vec4D m_l, m_n;
    inline LN_Pair(const ATOOLS::Vec4D &l=ATOOLS::Vec4D(),
		   const ATOOLS::Vec4D &n=ATOOLS::Vec4D(),
		   const int mode=0):
      m_l(mode&&l[0]<0.0?-l:l), m_n(mode&&n[0]<0.0?-n:n) {}
    inline ATOOLS::Vec4D SaneCMS() const
    { return (m_l[0]<0.0?-m_l:m_l)+(m_n[0]<0.0?-m_n:m_n); }
  };// end of struct LN_Pair

  LN_Pair GetLN(const ATOOLS::Vec4D &pi,
		const ATOOLS::Vec4D &pk,const int mode=0);

  ATOOLS::Vec4D LT(const ATOOLS::Vec4D &a,const ATOOLS::Vec4D &b,
		   const ATOOLS::Vec4D &c);

  double ComputePhi(ATOOLS::Vec4D pijt,ATOOLS::Vec4D pkt,
		    ATOOLS::Vec4D pi);

}// end of namespace PHASIC

#endif
