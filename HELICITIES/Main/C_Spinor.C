#include "C_Spinor.H"

#include "Exception.H"

#define ZERO Complex(0.0,0.0)
#define ONE Complex(1.0,0.0)

using namespace ATOOLS;
using namespace HELICITIES;

double CSpinor::s_accu(1.0e-12);

unsigned int CSpinor::s_r1(1), CSpinor::s_r2(2), CSpinor::s_r3(3);

std::ostream &HELICITIES::operator<<(std::ostream &ostr,const CSpinor &s)
{
  return ostr<<(s.B()>0?(s.R()>0?"|u,":"|v,"):(s.R()>0?"<u,":"<v,"))
	     <<s.H(0)<<","<<s.H(1)<<";"<<s()<<";"<<s[0]<<","<<s[1]<<","
	     <<s[2]<<","<<s[3]<<(s.B()>0?">":"|");
} 

std::ostream &HELICITIES::operator<<(std::ostream &ostr,
				 const std::vector<CSpinor> &v)
{
  ostr<<"{";
  for (int i(0);i<(int)v.size()-1;++i) ostr<<v[i]<<",";
  if (v.size()>0) ostr<<v.back();
  return ostr<<"}";
}

void CSpinor::SetGauge(const int gauge)
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

void CSpinor::Construct(const int h,const Vec4D &p,double m2)
{
  Vec4D ph(p[0]<0.0?-p.PSpat():p.PSpat(),p[1],p[2],p[3]);
  if (m_r>0^h<0) {// u+(p,m) / v-(p,m)
    Complex rpp(csqrt(PPlus(ph))), pt(PT(ph));
    m_u[2]=rpp;
    m_u[3]=IsZero(pt)?csqrt(PMinus(ph)):pt/rpp;
  }
  else {// u-(p,m) / v+(p,m)
    Complex rpp(csqrt(PPlus(ph))), pt(PTC(ph));
    m_u[0]=IsZero(pt)?csqrt(PMinus(ph)):pt/rpp;
    m_u[1]=-rpp;
  }
  if (m2<0.0) m2=p.Abs2();
  if (!IsZero(m2)) {
    double sgn(m_r>0?Sign(p[0]):-Sign(p[0]));
    double omp(sqrt(0.5*(dabs(p[0]/ph[0])+1.0)));
    double omm(sgn*sqrt(m2)/(2.0*dabs(ph[0])*omp));
    size_t r(m_r>0^h<0?0:2);
    m_u[0+r]=omm*m_u[2-r];
    m_u[1+r]=omm*m_u[3-r];
    m_u[2-r]*=omp;
    m_u[3-r]*=omp;
  }
  if (m_b<0) {
    m_b=1;
    *this=Bar();
  }
}

Complex CSpinor::operator*(const CSpinor &s) const
{ 
  if (s.m_b==m_b) THROW(fatal_error,"Equal spinor type");
  return m_u[0]*s.m_u[0]+m_u[1]*s.m_u[1]+m_u[2]*s.m_u[2]+m_u[3]*s.m_u[3];
}

bool CSpinor::operator==(const CSpinor &s) const
{
  double max(Max(std::abs(m_u[0]),
		 Max(std::abs(m_u[1]),
		     Max(std::abs(m_u[2]),std::abs(m_u[3]))))); 
  double q(IsZero(max)?1.0:1.0/max);
  for (short unsigned int i(0);i<4;++i) {
    if (std::abs(q*(m_u[i]-s.m_u[i]))>Accuracy()) return false;
  }
  return true;
}

CSpinor CSpinor::operator*(const double &d) const
{ 
  return CSpinor(m_r,m_b,m_u[0]*d,m_u[1]*d,m_u[2]*d,m_u[3]*d,
		 m_c,m_h[0],m_h[1]); 
}

CSpinor CSpinor::operator*(const Complex &c) const
{ 
  return CSpinor(m_r,m_b,m_u[0]*c,m_u[1]*c,m_u[2]*c,m_u[3]*c,
		 m_c,m_h[0],m_h[1]); 
}

CSpinor CSpinor::operator/(const double &d) const
{ 
  return CSpinor(m_r,m_b,m_u[0]/d,m_u[1]/d,m_u[2]/d,m_u[3]/d,
		 m_c,m_h[0],m_h[1]); 
}

CSpinor CSpinor::operator/(const Complex &c) const
{ 
  return CSpinor(m_r,m_b,m_u[0]/c,m_u[1]/c,m_u[2]/c,m_u[3]/c,
		 m_c,m_h[0],m_h[1]); 
}

CSpinor CSpinor::operator*=(const double &d) 
{ 
  m_u[0]*=d; 
  m_u[1]*=d; 
  m_u[2]*=d; 
  m_u[3]*=d; 
  return *this; 
}

CSpinor CSpinor::operator*=(const Complex &c) 
{ 
  m_u[0]*=c; 
  m_u[1]*=c; 
  m_u[2]*=c; 
  m_u[3]*=c; 
  return *this; 
}

CSpinor CSpinor::operator/=(const double &d) 
{ 
  m_u[0]/=d; 
  m_u[1]/=d; 
  m_u[2]/=d; 
  m_u[3]/=d; 
  return *this; 
}

CSpinor CSpinor::operator/=(const Complex &c) 
{ 
  m_u[0]/=c; 
  m_u[1]/=c; 
  m_u[2]/=c; 
  m_u[3]/=c; 
  return *this; 
}

CSpinor CSpinor::operator+(const CSpinor &s) const 
{ 
  return CSpinor(m_r,m_b,m_u[0]+s.m_u[0],m_u[1]+s.m_u[1],
		 m_u[2]+s.m_u[2],m_u[3]+s.m_u[3],m_c,m_h[0],m_h[1]); 
}

CSpinor CSpinor::operator-(const CSpinor &s) const
{ 
  return CSpinor(m_r,m_b,m_u[0]-s.m_u[0],m_u[1]-s.m_u[1],
		 m_u[2]-s.m_u[2],m_u[3]-s.m_u[3],m_c,m_h[0],m_h[1]); 
}

CSpinor CSpinor::operator+=(const CSpinor &s) 
{ 
  m_u[0]+=s.m_u[0]; 
  m_u[1]+=s.m_u[1]; 
  m_u[2]+=s.m_u[2]; 
  m_u[3]+=s.m_u[3]; 
  return *this; 
}

CSpinor CSpinor::operator-=(const CSpinor &s) 
{ 
  m_u[0]-=s.m_u[0]; 
  m_u[1]-=s.m_u[1]; 
  m_u[2]-=s.m_u[2]; 
  m_u[3]-=s.m_u[3]; 
  return *this; 
}

