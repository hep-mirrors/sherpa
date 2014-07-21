#include "METOOLS/Currents/C_Spinor.H"

#include "ATOOLS/Org/Exception.H"

#define ZERO SComplex(0.0,0.0)
#define ONE SComplex(1.0,0.0)

using namespace METOOLS;
using namespace ATOOLS;

template <class Scalar>
double CSpinor<Scalar>::s_accu(1.0e-12);

template <class Scalar> std::ostream &
METOOLS::operator<<(std::ostream &ostr,const CSpinor<Scalar> &s)
{
  return ostr<<(s.B()>0?(s.R()>=0?(s.R()==0?"|X(":"|u("):"|v("):
		(s.R()>=0?(s.R()==0?"<X(":"<u("):"<v("))<<s.On()
	     <<"),"<<s.H()<<","<<s.S()<<";"<<s(0)<<","<<s(1)<<";"
	     <<s[0]<<","<<s[1]<<","<<s[2]<<","<<s[3]<<(s.B()>0?">":"|");
} 

template <class Scalar> void CSpinor<Scalar>::
Construct(const int h,const Vec4<Scalar> &p,Scalar m2,const int ms)
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
    Spinor<Scalar> sh(1,ph); 
    m_u[2]=sh[0]; 
    m_u[3]=sh[1]; 
    m_on=2;
  } 
  else {// u-(p,m) / v+(p,m) 
    Spinor<Scalar> sh(-1,ph); 
    if (p[0]<0.0) sh=-sh;
    m_u[0]=sh[1]; 
    m_u[1]=-sh[0]; 
    m_on=1; 
  } 
  if (m2<0.0) m2=p.Abs2();
  if (!ATOOLS::IsZero(m2)) {
    Scalar sgn((m_r>0)^(ms<0)?Scalar(1.0):Scalar(-1.0));
    Scalar omp(sqrt((p[0]+ph[0])/(2.0*ph[0])));
    Scalar omm(sqrt((p[0]-ph[0])/(2.0*ph[0])));
    size_t r((m_r>0)^(h<0)?0:2);
    m_u[0+r]=sgn*omm*m_u[2-r];
    m_u[1+r]=sgn*omm*m_u[3-r];
    m_u[2-r]*=omp;
    m_u[3-r]*=omp;
    m_on=3;
  }
  if (abs(m_r)==2) m_r=0;
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
  if (s.m_b==m_b) { 
    if (m_r==0 && s.m_r==0) 
      return (*this)*s.CConj();
    THROW(fatal_error,"Equal spinor type");
  }
  return m_u[0]*s.m_u[0]+m_u[1]*s.m_u[1]+m_u[2]*s.m_u[2]+m_u[3]*s.m_u[3];
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::CConj() const
{ 
  if (m_b<0)
  return CSpinor(-m_r,-m_b,m_u[1],-m_u[0],-m_u[3],m_u[2],
		 m_c[0],m_c[1],m_h,m_s,m_on);
  return CSpinor(-m_r,-m_b,-m_u[1],m_u[0],m_u[3],-m_u[2],
		 m_c[0],m_c[1],m_h,m_s,m_on);
}

template <class Scalar>
bool CSpinor<Scalar>::operator==(const CSpinor<Scalar> &s) const
{
  Scalar max(Max(std::abs(m_u[0]),
		 Max(std::abs(m_u[1]),
		     Max(std::abs(m_u[2]),std::abs(m_u[3]))))); 
  Scalar q(ATOOLS::IsZero(max)?Scalar(1.0):Scalar(1.0)/max);
  for (short unsigned int i(0);i<4;++i) {
    if (ATOOLS::Abs(q*(m_u[i]-s.m_u[i]))>Accuracy()) return false;
  }
  return true;
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator*(const Scalar &d) const
{ 
  if (m_on==1)
  return CSpinor(m_r,m_b,m_u[0]*d,m_u[1]*d,ZERO,ZERO,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
  if (m_on==2)
  return CSpinor(m_r,m_b,ZERO,ZERO,m_u[2]*d,m_u[3]*d,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
  return CSpinor(m_r,m_b,m_u[0]*d,m_u[1]*d,m_u[2]*d,m_u[3]*d,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator*(const SComplex &c) const
{ 
  if (m_on==1)
  return CSpinor(m_r,m_b,m_u[0]*c,m_u[1]*c,ZERO,ZERO,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
  if (m_on==2)
  return CSpinor(m_r,m_b,ZERO,ZERO,m_u[2]*c,m_u[3]*c,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
  return CSpinor(m_r,m_b,m_u[0]*c,m_u[1]*c,m_u[2]*c,m_u[3]*c,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator/(const Scalar &d) const
{ 
  if (m_on==1)
  return CSpinor(m_r,m_b,m_u[0]/d,m_u[1]/d,ZERO,ZERO,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
  if (m_on==2)
  return CSpinor(m_r,m_b,ZERO,ZERO,m_u[2]/d,m_u[3]/d,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
  return CSpinor(m_r,m_b,m_u[0]/d,m_u[1]/d,m_u[2]/d,m_u[3]/d,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator/(const SComplex &c) const
{ 
  if (m_on==1)
  return CSpinor(m_r,m_b,m_u[0]/c,m_u[1]/c,ZERO,ZERO,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
  if (m_on==2)
  return CSpinor(m_r,m_b,ZERO,ZERO,m_u[2]/c,m_u[3]/c,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
  return CSpinor(m_r,m_b,m_u[0]/c,m_u[1]/c,m_u[2]/c,m_u[3]/c,
		 m_c[0],m_c[1],m_h,m_s,m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator*=(const Scalar &d) 
{ 
  if (m_on&1) {
  m_u[0]*=d; 
  m_u[1]*=d; 
  }
  if (m_on&2) {
  m_u[2]*=d; 
  m_u[3]*=d; 
  } 
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator*=(const SComplex &c) 
{ 
  if (m_on&1) {
  m_u[0]*=c; 
  m_u[1]*=c; 
  }
  if (m_on&2) {
  m_u[2]*=c; 
  m_u[3]*=c;
  } 
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator/=(const Scalar &d) 
{ 
  if (m_on&1) {
  m_u[0]/=d; 
  m_u[1]/=d; 
  }
  if (m_on&2) {
  m_u[2]/=d; 
  m_u[3]/=d; 
  } 
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator/=(const SComplex &c) 
{ 
  if (m_on&1) {
  m_u[0]/=c; 
  m_u[1]/=c; 
  }
  if (m_on&2) {
  m_u[2]/=c; 
  m_u[3]/=c; 
  } 
  return *this; 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator+(const CSpinor &s) const 
{ 
  return CSpinor(m_r,m_b,m_u[0]+s.m_u[0],m_u[1]+s.m_u[1],
		 m_u[2]+s.m_u[2],m_u[3]+s.m_u[3],m_c[0],m_c[1],
		 m_h,m_s,m_on|s.m_on); 
}

template <class Scalar>
CSpinor<Scalar> CSpinor<Scalar>::operator-(const CSpinor &s) const
{ 
  return CSpinor(m_r,m_b,m_u[0]-s.m_u[0],m_u[1]-s.m_u[1],
		 m_u[2]-s.m_u[2],m_u[3]-s.m_u[3],m_c[0],m_c[1],
		 m_h,m_s,m_on|s.m_on); 
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

template <class Scalar>
void CSpinor<Scalar>::Add(const CObject *c)
{
  const CSpinor *s(static_cast<const CSpinor*>(c));
  m_on|=s->m_on;
  if (m_on&1) {
  m_u[0]+=s->m_u[0]; 
  m_u[1]+=s->m_u[1]; 
  }
  if (m_on&2) {
  m_u[2]+=s->m_u[2]; 
  m_u[3]+=s->m_u[3]; 
  }
}

template <class Scalar>
void CSpinor<Scalar>::Divide(const double &d)
{
  m_u[0]/=Scalar(d);
  m_u[1]/=Scalar(d); 
  m_u[2]/=Scalar(d); 
  m_u[3]/=Scalar(d);
}

template <class Scalar>
void CSpinor<Scalar>::Invert()
{
  m_u[0]=-m_u[0];
  m_u[1]=-m_u[1];
  m_u[2]=-m_u[2];
  m_u[3]=-m_u[3];
}

template <class Scalar>
bool CSpinor<Scalar>::IsZero() const
{
  return m_u[0]==Scalar(0.0) &&
    m_u[1]==Scalar(0.0) &&
    m_u[2]==Scalar(0.0) &&
    m_u[3]==Scalar(0.0);
}

template <class Scalar>
typename ATOOLS::AutoDelete_Vector<CSpinor<Scalar> > 
CSpinor<Scalar>::s_objects;

template <class Scalar>
CSpinor<Scalar> *CSpinor<Scalar>::New()
{
  if (s_objects.empty()) {
    return new CSpinor();
  }
  CSpinor *v(s_objects.back());
  s_objects.pop_back();
  return v;
}

template <class Scalar>
CSpinor<Scalar> *CSpinor<Scalar>::New(const CSpinor &s)
{
  if (s_objects.empty()) {
    return new CSpinor(s);
  }
  CSpinor *v(s_objects.back());
  s_objects.pop_back();
  *v=s;
  return v;
}

template <class Scalar>
CSpinor<Scalar> *CSpinor<Scalar>::New
(const int r,const int b,const int cr,const int ca,
 const size_t &h,const size_t &s,const int on)
{
  if (s_objects.empty()) {
    return new CSpinor(r,b,cr,ca,h,s,on);
  }
  CSpinor *v(s_objects.back());
  s_objects.pop_back();
  v->m_r=r;
  v->m_b=b;
  v->m_on=on;
  v->m_u[0]=v->m_u[1]=v->m_u[2]=v->m_u[3]=Scalar(0.0);
  v->m_h=h;
  v->m_s=s;
  v->m_c[0]=cr;
  v->m_c[1]=ca;
  return v;
}

template <class Scalar>
CObject *CSpinor<Scalar>::Copy() const
{
  return CSpinor::New(*this);
}

template <class Scalar>
void CSpinor<Scalar>::Delete()
{
  s_objects.push_back(this);
}

namespace METOOLS {

  template class DDSpinor;
  template std::ostream &operator<<(std::ostream &ostr,const DDSpinor &s);

  template class QDSpinor;
  template std::ostream &operator<<(std::ostream &ostr,const QDSpinor &s);

}
