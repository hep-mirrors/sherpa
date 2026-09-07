#include "ATOOLS/Phys/Spinor.H"

#include "ATOOLS/Org/Exception.H"

#include <cstdlib>

using namespace ATOOLS;

// #define TEST_Representation

namespace {
  // The threshold below which Spinor::Construct() abandons the accurate
  // m_u2 = pT/sqrt(p+) form and falls back to sqrt(p-). That fallback is a
  // hard branch, not a rounding effect, so it survives any increase in
  // precision - which makes it the first thing to scan when higher precision
  // does not move the answer. SHERPA_SPINOR_ACCU overrides it for that scan.
  double SpinorAccuDefault()
  {
    const char *e(getenv("SHERPA_SPINOR_ACCU"));
    return e?atof(e):1.0e-12;
  }
}

template <class Scalar>
double Spinor<Scalar>::s_accu(SpinorAccuDefault());

template <class Scalar> std::ostream &
ATOOLS::operator<<(std::ostream &ostr,const Spinor<Scalar> &s)
{
  return ostr<<"|"<<s(0)<<","<<s(1)<<(s.R()>0?">":"]");
} 

template <class Scalar>
void Spinor<Scalar>::SetGauge(const int gauge)
{
  switch (gauge) {
  case 0: s_r1=1; s_r2=2; s_r3=3; break;
  case 1: s_r1=2; s_r2=3; s_r3=1; break;
  case 2: s_r1=3; s_r2=1; s_r3=2; break;
  default:
    THROW(fatal_error,"Gauge choice not implemented");
  }
}

template <class Scalar> Vec4<Scalar> Spinor<Scalar>::GetK0()
{
  Vec4<Scalar> k0(1.0,0.0,0.0,0.0);
  k0[R3()]=-1.0;
  return k0;
}

template <class Scalar> Vec4<Scalar> Spinor<Scalar>::GetK1()
{
  Vec4<Scalar> k1(0.0,0.0,0.0,0.0);
  k1[R1()]=1.0;
  return k1;
}

template <class Scalar>
void Spinor<Scalar>::Construct(const Vec4<Scalar> &p)
{
  // Scalar/SComplex throughout, not double/Complex: the light-cone
  // components and their square roots are the whole content of the spinor,
  // and rounding them to double here would discard every extra digit a wider
  // Scalar was instantiated to provide.
  Scalar pp(PPlus(p)), pm(PMinus(p));
  SComplex rpp(csqrt(pp)), rpm(csqrt(pm)), pt(PT(p));
  m_u1=rpp;
  m_u2=rpm;
  Scalar sv(Abs(p[0])*Scalar(s_accu));
  if ((Abs(pt.real())>sv || Abs(pt.imag())>sv) &&
      (Abs(rpp.real())>sv || Abs(rpp.imag())>sv)) {
    m_u2=SComplex(pt.real(),m_r>0?pt.imag():-pt.imag())/rpp;
  }
#ifndef USING__old_phase_convention
  if (pp<Scalar(0.0) || pm<Scalar(0.0)) {
    if (m_r<0) {
      m_u1=SComplex(-m_u1.imag(),m_u1.real());
      m_u2=SComplex(-m_u2.imag(),m_u2.real());
    }
    else {
      m_u1=-SComplex(-m_u1.imag(),m_u1.real());
      m_u2=-SComplex(-m_u2.imag(),m_u2.real());
    }
  }
#endif
}

template <class Scalar>
void Spinor<Scalar>::ConstructLC(const Vec4<Scalar> &p)
{
  // Same construction as Construct(), except for how the two light-cone
  // components are obtained. Directly, one of p[0]+p[3] and p[0]-p[3] is a
  // difference of two nearly equal numbers: for a momentum collinear to the
  // beam it loses every digit (measured at theta = 5e-5, w = 1e-5 GeV it
  // returns exactly 0.0 where the true value is 2.7e-21), and the spinor is
  // then built from a zero. Since p^2 = 0 here, (p0+pz)(p0-pz) = pT^2 holds
  // exactly, so the small component follows from the large one by a division
  // with no cancellation anywhere. Accurate to ~6e-17 relative at the same
  // kinematics.
  Scalar pp(PPlus(p)), pm(PMinus(p));
  const Scalar pt2(p[s_r1]*p[s_r1]+p[s_r2]*p[s_r2]);
  // Divide by whichever has the LARGER MAGNITUDE: pp*pm = pT^2 >= 0 so the two
  // always share a sign, and for p[0]<0 (CSpinor passes -PSpat() there) both
  // are negative -- comparing them with > rather than by magnitude would pick
  // the cancelling one and divide by the very quantity that lost its digits.
  if (dabs(pp)>dabs(pm)) { if (pp!=Scalar(0.0)) pm=pt2/pp; }
  else                   { if (pm!=Scalar(0.0)) pp=pt2/pm; }
  SComplex rpp(csqrt(pp)), rpm(csqrt(pm)), pt(PT(p));
  m_u1=rpp;
  m_u2=rpm;
  Scalar sv(Abs(p[0])*Scalar(s_accu));
  if ((Abs(pt.real())>sv || Abs(pt.imag())>sv) &&
      (Abs(rpp.real())>sv || Abs(rpp.imag())>sv)) {
    m_u2=SComplex(pt.real(),m_r>0?pt.imag():-pt.imag())/rpp;
  }
#ifndef USING__old_phase_convention
  if (pp<Scalar(0.0) || pm<Scalar(0.0)) {
    if (m_r<0) {
      m_u1=SComplex(-m_u1.imag(),m_u1.real());
      m_u2=SComplex(-m_u2.imag(),m_u2.real());
    }
    else {
      m_u1=-SComplex(-m_u1.imag(),m_u1.real());
      m_u2=-SComplex(-m_u2.imag(),m_u2.real());
    }
  }
#endif
}

template <class Scalar> 
std::complex<Scalar> Spinor<Scalar>::operator*(const Spinor &s) const
{ 
#ifdef TEST_Representation
  if (m_r!=s.m_r) {
    msg_Error()<<METHOD<<"(..): Distinct representations."<<std::endl;
    return Complex(0.0,0.0);
  }
#endif
  return m_u1*s.m_u2-m_u2*s.m_u1; 
}

template <class Scalar>
bool Spinor<Scalar>::operator==(const Spinor &s) const
{
  Scalar max(Max(Abs(m_u1),Abs(m_u2)));
  Scalar q(IsZero(max)?Scalar(1.0):Scalar(1.0)/max);
  if (Abs(q*(m_u1-s.m_u1))>Scalar(Accuracy())) return false;
  if (Abs(q*(m_u2-s.m_u2))>Scalar(Accuracy())) return false;
  return true;
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator*(const Scalar &d) const
{ 
  return Spinor(m_r,m_u1*d,m_u2*d); 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator*(const SComplex &c) const
{ 
  return Spinor(m_r,m_u1*c,m_u2*c); 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator/(const Scalar &d) const
{ 
  return Spinor(m_r,m_u1/d,m_u2/d); 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator/(const SComplex &c) const
{ 
  return Spinor(m_r,m_u1/c,m_u2/c); 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator*=(const Scalar &d) 
{ 
  m_u1*=d; 
  m_u2*=d; 
  return *this; 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator*=(const SComplex &c) 
{ 
  m_u1*=c; 
  m_u2*=c; 
  return *this; 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator/=(const Scalar &d) 
{ 
  m_u1/=d; 
  m_u2/=d; 
  return *this; 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator/=(const SComplex &c) 
{ 
  m_u1/=c;
  m_u2/=c; 
  return *this; 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator+(const Spinor &s) const 
{ 
  return Spinor(m_r,m_u1+s.m_u1,m_u2+s.m_u2); 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator-(const Spinor &s) const
{ 
  return Spinor(m_r,m_u1-s.m_u1,m_u2-s.m_u2); 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator+=(const Spinor &s) 
{ 
  m_u1+=s.m_u1; 
  m_u2+=s.m_u2; 
  return *this; 
}

template <class Scalar>
Spinor<Scalar> Spinor<Scalar>::operator-=(const Spinor &s) 
{ 
  m_u1-=s.m_u1; 
  m_u2-=s.m_u2; 
  return *this; 
}

namespace ATOOLS {

  template class DWSpinor;
  template std::ostream &operator<<(std::ostream &ostr,const DWSpinor &s);

  template class QWSpinor;
  template std::ostream &operator<<(std::ostream &ostr,const QWSpinor &s);

  template class XWSpinor;
  template std::ostream &operator<<(std::ostream &ostr,const XWSpinor &s);

}
