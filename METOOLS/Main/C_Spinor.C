#include "METOOLS/Main/C_Spinor.H"

#include "ATOOLS/Org/Exception.H"

#define ZERO SComplex(0.0,0.0)
#define ONE SComplex(1.0,0.0)

using namespace ATOOLS;

template <class Scalar>
double CSpinor<Scalar>::s_accu(1.0e-12);

template <class Scalar>
unsigned int CSpinor<Scalar>::s_r1(1);
template <class Scalar>
unsigned int CSpinor<Scalar>::s_r2(2);
template <class Scalar>
unsigned int CSpinor<Scalar>::s_r3(3);

template <class Scalar> std::ostream &
ATOOLS::operator<<(std::ostream &ostr,const CSpinor<Scalar> &s)
{
  return ostr<<(s.B()>0?(s.R()>0?"|u(":"|v("):(s.R()>0?"<u(":"<v("))<<s.On()
	     <<"),"<<s.H(0)<<","<<s.H(1)<<";"<<s()<<";"<<s[0]<<","<<s[1]<<","
	     <<s[2]<<","<<s[3]<<(s.B()>0?">":"|");
} 

template <class Scalar> std::ostream &
ATOOLS::operator<<(std::ostream &ostr,
		   const std::vector<CSpinor<Scalar> > &v)
{
  ostr<<"{";
  for (int i(0);i<(int)v.size()-1;++i) ostr<<v[i]<<",";
  if (v.size()>0) ostr<<v.back();
  return ostr<<"}";
}

template <class Scalar>
void CSpinor<Scalar>::SetGauge(const int gauge)
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

template <class Scalar> void CSpinor<Scalar>::
Construct(const int h,const Vec4<Scalar> &p,Scalar m2)
{
  if (Abs(p[1])==0.0 && Abs(p[2])==0.0 && Abs(p[3])==0.0) {
    SComplex rte(csqrt(p[0]));
    if ((m_r>0)^(h<0)) {// u+(p,m) / v-(p,m)
      m_u[2]=rte;
      m_u[3]=ZERO;
    }
    else {// u-(p,m) / v+(p,m)
      m_u[0]=ZERO;
      m_u[1]=-rte;
    }
    Scalar sgn(m_r>0?Scalar(1.0):Scalar(-1.0));
    size_t r((m_r>0)^(h<0)?0:2);
    m_u[0+r]=sgn*m_u[2-r];
    m_u[1+r]=sgn*m_u[3-r];
    m_on=3;
    if (m_b<0) {
      m_b=1;
      *this=Bar();
    }
  }
  else {
  Vec4<Scalar> ph(p[0]<0.0?-p.PSpat():p.PSpat(),p[1],p[2],p[3]);
  if ((m_r>0)^(h<0)) {// u+(p,m) / v-(p,m)
    SComplex rpp(csqrt(PPlus(ph))), pt(PT(ph));
    m_u[2]=rpp;
    m_u[3]=IsZero(pt)?csqrt(PMinus(ph)):pt/rpp;
    m_on=2;
  }
  else {// u-(p,m) / v+(p,m)
    SComplex rpp(csqrt(PPlus(ph))), pt(PTC(ph));
    m_u[0]=IsZero(pt)?csqrt(PMinus(ph)):pt/rpp;
    m_u[1]=-rpp;
    m_on=1;
  }
  if (m2<0.0) m2=p.Abs2();
  if (!IsZero(m2)) {
    Scalar sgn(m_r>0?Scalar(1.0):Scalar(-1.0));
    Scalar omp(sqrt((p[0]+ph[0])/(2.0*ph[0])));
    Scalar omm(sqrt((p[0]-ph[0])/(2.0*ph[0])));
    size_t r((m_r>0)^(h<0)?0:2);
    m_u[0+r]=sgn*omm*m_u[2-r];
    m_u[1+r]=sgn*omm*m_u[3-r];
    m_u[2-r]*=omp;
    m_u[3-r]*=omp;
    m_on=3;
  }
  if (m_b<0) {
    m_b=1;
    *this=Bar();
  }
  }
}

template <class Scalar> bool CSpinor<Scalar>::SetOn()
{
  m_on=0;
  if (m_u[0]!=ZERO || m_u[1]!=ZERO) m_on|=1;
  if (m_u[2]!=ZERO || m_u[3]!=ZERO) m_on|=2;
  return m_on&3;
}

template <class Scalar> std::complex<Scalar> 
CSpinor<Scalar>::operator*(const CSpinor<Scalar> &s) const
{ 
  if (s.m_b==m_b) THROW(fatal_error,"Equal spinor type");
  return m_u[0]*s.m_u[0]+m_u[1]*s.m_u[1]+m_u[2]*s.m_u[2]+m_u[3]*s.m_u[3];
}

template <class Scalar>
bool CSpinor<Scalar>::operator==(const CSpinor<Scalar> &s) const
{
  Scalar max(Max(std::abs(m_u[0]),
		 Max(std::abs(m_u[1]),
		     Max(std::abs(m_u[2]),std::abs(m_u[3]))))); 
  Scalar q(IsZero(max)?Scalar(1.0):Scalar(1.0)/max);
  for (short unsigned int i(0);i<4;++i) {
    if (ATOOLS::Abs(q*(m_u[i]-s.m_u[i]))>Accuracy()) return false;
  }
  return true;
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator*(const Scalar &d) const
{ 
  return CSpinor(m_r,m_b,m_u[0]*d,m_u[1]*d,m_u[2]*d,m_u[3]*d,
		 m_c,m_h[0],m_h[1],m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator*(const SComplex &c) const
{ 
  return CSpinor(m_r,m_b,m_u[0]*c,m_u[1]*c,m_u[2]*c,m_u[3]*c,
		 m_c,m_h[0],m_h[1],m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator/(const Scalar &d) const
{ 
  return CSpinor(m_r,m_b,m_u[0]/d,m_u[1]/d,m_u[2]/d,m_u[3]/d,
		 m_c,m_h[0],m_h[1],m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator/(const SComplex &c) const
{ 
  return CSpinor(m_r,m_b,m_u[0]/c,m_u[1]/c,m_u[2]/c,m_u[3]/c,
		 m_c,m_h[0],m_h[1],m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator*=(const Scalar &d) 
{ 
  m_u[0]*=d; 
  m_u[1]*=d; 
  m_u[2]*=d; 
  m_u[3]*=d; 
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator*=(const SComplex &c) 
{ 
  m_u[0]*=c; 
  m_u[1]*=c; 
  m_u[2]*=c; 
  m_u[3]*=c; 
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator/=(const Scalar &d) 
{ 
  m_u[0]/=d; 
  m_u[1]/=d; 
  m_u[2]/=d; 
  m_u[3]/=d; 
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator/=(const SComplex &c) 
{ 
  m_u[0]/=c; 
  m_u[1]/=c; 
  m_u[2]/=c; 
  m_u[3]/=c; 
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator+(const CSpinor &s) const 
{ 
  return CSpinor(m_r,m_b,m_u[0]+s.m_u[0],m_u[1]+s.m_u[1],
		 m_u[2]+s.m_u[2],m_u[3]+s.m_u[3],m_c,m_h[0],m_h[1],
		 m_on|s.m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator-(const CSpinor &s) const
{ 
  return CSpinor(m_r,m_b,m_u[0]-s.m_u[0],m_u[1]-s.m_u[1],
		 m_u[2]-s.m_u[2],m_u[3]-s.m_u[3],m_c,m_h[0],m_h[1],
		 m_on|s.m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator+=(const CSpinor &s) 
{ 
  m_u[0]+=s.m_u[0]; 
  m_u[1]+=s.m_u[1]; 
  m_u[2]+=s.m_u[2]; 
  m_u[3]+=s.m_u[3]; 
  m_on|=s.m_on;
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator-=(const CSpinor &s) 
{ 
  m_u[0]-=s.m_u[0]; 
  m_u[1]-=s.m_u[1]; 
  m_u[2]-=s.m_u[2]; 
  m_u[3]-=s.m_u[3]; 
  m_on|=s.m_on;
  return *this; 
}

namespace ATOOLS {

  template class DDSpinor;
  template std::ostream &operator<<(std::ostream &ostr,const DDSpinor &s);

  template class QDSpinor;
  template std::ostream &operator<<(std::ostream &ostr,const QDSpinor &s);

}
