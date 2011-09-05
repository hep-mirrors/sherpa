#include "ATOOLS/Phys/Spinor.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

// #define TEST_Representation

template <class Scalar>
double Spinor<Scalar>::s_accu(1.0e-12);

template <class Scalar>
unsigned int Spinor<Scalar>::s_r1(1);
template <class Scalar>
unsigned int Spinor<Scalar>::s_r2(2);
template <class Scalar>
unsigned int Spinor<Scalar>::s_r3(3);

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
  case 3: s_r1=1; s_r2=3; s_r3=2; break;
  case 4: s_r1=3; s_r2=2; s_r3=1; break;
  case 5: s_r1=2; s_r2=1; s_r3=3; break;
  default:
    THROW(fatal_error,"Gauge choice not implemented");
  }
}

template <class Scalar>
void Spinor<Scalar>::Construct(const Vec4<Scalar> &p)
{
  Complex rpp(csqrt(PPlus(p))), rpm(csqrt(PMinus(p))), pt(PT(p));
  double accu(sqrt(rpa->gen.Accu()));
  if (((rpp==Complex(0.0,0.0) || rpm==Complex(0.0,0.0)) &&
       pt!=Complex(0.0,0.0)) || dabs(std::abs(pt/(rpp*rpm))-1.0)>accu) {
    msg_Error()<<METHOD<<"(): Warning: \\sqrt{p^+p^-} = "<<std::abs(rpp*rpm)
	       <<" vs. |p_\\perp| = "<<std::abs(pt)<<", rel. diff. "
	       <<(std::abs(pt/(rpp*rpm))-1.0)<<"."<<std::endl;
  }
  m_u1=rpp;
  m_u2=rpm;
  if (pt!=Complex(0.0,0.0)) { 
    m_u2=Complex(pt.real(),m_r>0?pt.imag():-pt.imag())/rpp;
  }
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
  Scalar max(Max(std::abs(m_u1),std::abs(m_u2))); 
  Scalar q(IsZero(max)?1.0:1.0/max);
  if (std::abs(q*(m_u1-s.m_u1))>Accuracy()) return false;
  if (std::abs(q*(m_u2-s.m_u2))>Accuracy()) return false;
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

}
