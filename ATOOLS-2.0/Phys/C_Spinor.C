#include "C_Spinor.H"

#include "Exception.H"

using namespace ATOOLS;

const double sqrttwo(sqrt(2.0));

double CSpinor::s_accu(1.0e-12);

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const CSpinor &s)
{
  return ostr<<(s.R()>0?"|":"<")<<s()<<","<<s[0]<<","<<s[1]<<","
	     <<s[2]<<","<<s[3]<<(s.R()>0?">":"|");
} 

std::ostream &ATOOLS::operator<<(std::ostream &ostr,
				 const std::vector<CSpinor> &v)
{
  ostr<<"{";
  for (int i(0);i<(int)v.size()-1;++i) ostr<<v[i]<<",";
  if (v.size()>0) ostr<<v.back();
  return ostr<<"}";
}

void CSpinor::Construct(const int h,const Vec4D &p)
{
  Complex rpp(csqrt(PPlus(p))), rpm(csqrt(PMinus(p))), pt(PT(p));
  static double accu(sqrt(Accu()));
  if (((rpp==Complex(0.0,0.0) || rpm==Complex(0.0,0.0)) &&
       pt!=Complex(0.0,0.0)) || dabs(std::abs(pt/(rpp*rpm))-1.0)>accu) {
    msg.Error()<<METHOD<<"(): \\sqrt{p^+p^-} = "<<std::abs(rpp*rpm)
	       <<" vs. |p_\\perp| = "<<std::abs(pt)<<", rel. diff. "
	       <<(std::abs(pt/(rpp*rpm))-1.0)<<"."<<std::endl;
    THROW(fatal_error,"Cannot construct massive two-component spinor.");
  }
  Complex u1(rpp/sqrttwo), u2(rpm/sqrttwo);
  if (pt!=Complex(0.0,0.0)) { 
    u2*=Complex(pt.real(),h>0?pt.imag():-pt.imag())/std::abs(pt);
  }
  if (m_r>0) {
    // v(p,h)
    if (h>0) {
      // v_+(p)
      m_u[2]=-(m_u[0]=u2);
      m_u[3]=-(m_u[1]=-u1);
    }
    else {
      // v_-(p)
      m_u[2]=m_u[0]=u1;
      m_u[3]=m_u[1]=u2;
    }
  }
  else {
    // \bar u(p,h)
    if (h>0) {
      // \bar u_+(p)
      m_u[2]=-(m_u[0]=u1);
      m_u[3]=-(m_u[1]=u2);
    }
    else {
      // \bar u_-(p)
      m_u[2]=m_u[0]=u2;
      m_u[3]=m_u[1]=-u1;
    }
  }
}

Complex CSpinor::operator*(const CSpinor &s) const
{ 
  if (s.m_r==m_r) THROW(fatal_error,"Equal spinor type");
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
  return CSpinor(m_r,m_u[0]*d,m_u[1]*d,m_u[2]*d,m_u[3]*d,m_c); 
}

CSpinor CSpinor::operator*(const Complex &c) const
{ 
  return CSpinor(m_r,m_u[0]*c,m_u[1]*c,m_u[2]*c,m_u[3]*c,m_c); 
}

CSpinor CSpinor::operator/(const double &d) const
{ 
  return CSpinor(m_r,m_u[0]/d,m_u[1]/d,m_u[2]/d,m_u[3]/d,m_c); 
}

CSpinor CSpinor::operator/(const Complex &c) const
{ 
  return CSpinor(m_r,m_u[0]/c,m_u[1]/c,m_u[2]/c,m_u[3]/c,m_c); 
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
  return CSpinor(m_r,m_u[0]+s.m_u[0],m_u[1]+s.m_u[1],
		 m_u[2]+s.m_u[2],m_u[3]+s.m_u[3],m_c); 
}

CSpinor CSpinor::operator-(const CSpinor &s) const
{ 
  return CSpinor(m_r,m_u[0]-s.m_u[0],m_u[1]-s.m_u[1],
		 m_u[2]-s.m_u[2],m_u[3]-s.m_u[3],m_c); 
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

