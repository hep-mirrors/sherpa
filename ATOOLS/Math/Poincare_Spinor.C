#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Poincare.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Main/Polarization_Tools.H"   

using namespace ATOOLS;


/*void BoostBack_Spinor(METOOLS::CSpinor<double> &s, const ATOOLS::Vec4D &l) 
{
  //TODO: Spinor boost back (spinor is not a 4vector - another representation of Lorentz group)
  // need rapidity
  double m(sqrt(l.Abs2()));
  if (m < 1e-12) return;  
  //if (s.B()<0) {
    //METOOLS::CSpinor<double> k(s.Bar());
    //BoostBack(k);
    //s=k.Bar();
    //return;
  //}
  double n(1.0/sqrt(2.0*m*(l[0]+m)));
  // cosh(eta)
  double c(n*(l[0]+m));
  // sinh * vec(n)
  double sx(n*l[1]), sy(n*l[2]), sz(n*l[3]);
  double sgn(-1.0);
  for (int i(0);i<4;i+=2,sgn=-sgn) {
    Complex a(s[i]), b(s[i+1]);
    s[i] = c*a+sgn*(sz*a+Complex(sx,-sy)*b);
    s[i+1] = c*b+sgn*(Complex(sx,sy)*a-sz*b);
  }

}

void METOOLS::Polarization_Spinor_Basis::BoostBack_s(METOOLS::CSpinor<double>& s) {
  ATOOLS::Poincare_Spinor spinboost(beta);
  for (size_t i(0); i<m_s.size(); ++i) spinboost.BoostBack(m_s[i]);
}*/