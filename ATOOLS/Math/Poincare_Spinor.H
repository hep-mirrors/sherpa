#ifndef ATOOLS_Math_Poincare_Spinor_H
#define ATOOLS_Math_Poincare_Spinor_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "ATOOLS/Phys/Spinor.H"
#include <type_traits>

namespace ATOOLS {

  class Spinor_Boost {
  private:
  
    //int m_type;
    //ATOOLS::Vec4D  m_l;
    double m_c, m_sx, m_sy, m_sz;
    template<class S> void Apply(S &s, double sgn0) const{
      double sgn(sgn0);
      for (int i(0); i<4; i+=2, sgn=-sgn) {
        Complex a(s[i]), b(s[i+1]);
        s[i] = m_c*a + sgn*(m_sz*a + Complex(m_sx, -m_sy)*b);
        s[i+1] = m_c*b + sgn*(Complex(m_sx,m_sy)*a - m_sz*b);
      }
    }

  public:
    Spinor_Boost(const Vec4D &l) {
      const double m2(l.Abs2());
      //if (!(m2 > 1.0e-12) || !(l[0] > 0.0)) {
        //THROW(fatal_error, "Spinor_Boost needs a time-like four-momentum with E>0")
      //}
      const double m(std::sqrt(m2));
      const double n(1.0/std::sqrt(2.0*m*(l[0]+m)));
      m_c  = n*(l[0]+m);
      m_sx = n*l[1]; m_sy = n*l[2]; m_sz = n*l[3];
    }

    // Ruhesystem von l -> Lab (wie Poincare::BoostBack)
    template <class S> void BoostBack(S &s) const { Apply(s, -1.0); }
    // Lab -> Ruhesystem von l (Inverse, wie Poincare::Boost)
    template <class S> void Boost(S &s) const { Apply(s, +1.0); }
  };
}// end of namespace ATOOLS

#endif
