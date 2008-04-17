#include "Spinor.H"

#include "Exception.H"

using namespace ATOOLS;

// #define TEST_Representation

double Spinor::s_accu(1.0e-12);

unsigned int Spinor::s_r1(1), Spinor::s_r2(2), Spinor::s_r3(3);

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const Spinor &s)
{
  return ostr<<"|"<<s(0)<<","<<s(1)<<(s.R()>0?">":"]");
} 

std::ostream &ATOOLS::operator<<(std::ostream &ostr,
			      const std::vector<Spinor> &v)
{
  ostr<<"{";
  for (int i(0);i<(int)v.size()-1;++i) ostr<<v[i]<<",";
  if (v.size()>0) ostr<<v.back();
  return ostr<<"}";
}

void Spinor::SetGauge(const int gauge)
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

void Spinor::Construct(const ATOOLS::Vec4D &p)
{
  Complex rpp(csqrt(PPlus(p))), rpm(csqrt(PMinus(p))), pt(PT(p));
  static double accu(sqrt(Accu()));
  if (((rpp==Complex(0.0,0.0) || rpm==Complex(0.0,0.0)) &&
       pt!=Complex(0.0,0.0)) || dabs(std::abs(pt/(rpp*rpm))-1.0)>accu) {
    msg_Error()<<METHOD<<"(): \\sqrt{p^+p^-} = "<<std::abs(rpp*rpm)
	       <<" vs. |p_\\perp| = "<<std::abs(pt)<<", rel. diff. "
	       <<(std::abs(pt/(rpp*rpm))-1.0)<<"."<<std::endl;
    THROW(fatal_error,"Cannot construct massive two-component spinor.");
  }
  m_u1=rpp;
  m_u2=rpm;
  if (pt!=Complex(0.0,0.0)) { 
    m_u2=Complex(pt.real(),m_r>0?pt.imag():-pt.imag())/rpp;
  }
}

Complex Spinor::operator*(const Spinor &s) const
{ 
#ifdef TEST_Representation
  if (m_r!=s.m_r) {
    msg_Error()<<METHOD<<"(..): Distinct representations."<<std::endl;
    return Complex(0.0,0.0);
  }
#endif
  return m_u1*s.m_u2-m_u2*s.m_u1; 
}

bool Spinor::operator==(const Spinor &s) const
{
  double max(Max(std::abs(m_u1),std::abs(m_u2))); 
  double q(IsZero(max)?1.0:1.0/max);
  if (std::abs(q*(m_u1-s.m_u1))>Accuracy()) return false;
  if (std::abs(q*(m_u2-s.m_u2))>Accuracy()) return false;
  return true;
}

Spinor Spinor::operator*(const double &d) const
{ 
  return Spinor(m_r,m_u1*d,m_u2*d); 
}

Spinor Spinor::operator*(const Complex &c) const
{ 
  return Spinor(m_r,m_u1*c,m_u2*c); 
}

Spinor Spinor::operator/(const double &d) const
{ 
  return Spinor(m_r,m_u1/d,m_u2/d); 
}

Spinor Spinor::operator/(const Complex &c) const
{ 
  return Spinor(m_r,m_u1/c,m_u2/c); 
}

Spinor Spinor::operator*=(const double &d) 
{ 
  m_u1*=d; 
  m_u2*=d; 
  return *this; 
}

Spinor Spinor::operator*=(const Complex &c) 
{ 
  m_u1*=c; 
  m_u2*=c; 
  return *this; 
}

Spinor Spinor::operator/=(const double &d) 
{ 
  m_u1/=d; 
  m_u2/=d; 
  return *this; 
}

Spinor Spinor::operator/=(const Complex &c) 
{ 
  m_u1/=c;
  m_u2/=c; 
  return *this; 
}

Spinor Spinor::operator+(const Spinor &s) const 
{ 
  return Spinor(m_r,m_u1+s.m_u1,m_u2+s.m_u2); 
}

Spinor Spinor::operator-(const Spinor &s) const
{ 
  return Spinor(m_r,m_u1-s.m_u1,m_u2-s.m_u2); 
}

Spinor Spinor::operator+=(const Spinor &s) 
{ 
  m_u1+=s.m_u1; 
  m_u2+=s.m_u2; 
  return *this; 
}

Spinor Spinor::operator-=(const Spinor &s) 
{ 
  m_u1-=s.m_u1; 
  m_u2-=s.m_u2; 
  return *this; 
}

