#include "METOOLS/Currents/C_Pseudo.H"

using namespace METOOLS;
using namespace ATOOLS;

template<class Scalar>
double CAsT4<Scalar>::s_accu(1.0e-12);

template<class Scalar> std::ostream &
METOOLS::operator<<(std::ostream &s,const CAsT4<Scalar> &ten)
{
  return s<<'{'<<ten.H()<<","<<ten.S()<<";"<<ten(0)<<","<<ten(1)<<'|'
	  <<ten[0]<<','<<ten[1]<<','<<ten[2]<<','
	  <<ten[3]<<','<<ten[4]<<','<<ten[5]<<'}';
}

template<class Scalar>
bool CAsT4<Scalar>::Nan() const
{
  for(short unsigned int i(0);i<6;++i) {
    if (IsNan(m_x[i])) return true;
  }
  return false;
}

template<class Scalar>
void CAsT4<Scalar>::ResetAccu()                
{ 
  s_accu=Accu(); 
}

template<class Scalar>
CAsT4<Scalar>::CAsT4(const CVec4<Scalar> &v1,const CVec4<Scalar> &v2)
{ 
  m_x[0]=v1[0]*v2[1]-v1[1]*v2[0];// 0,1 
  m_x[1]=v1[1]*v2[2]-v1[2]*v2[1];// 1,2 
  m_x[2]=v1[2]*v2[3]-v1[3]*v2[2];// 2,3 
  m_x[3]=v1[0]*v2[2]-v1[2]*v2[0];// 0,2 
  m_x[4]=v1[1]*v2[3]-v1[3]*v2[1];// 1,3 
  m_x[5]=v1[0]*v2[3]-v1[3]*v2[0];// 0,3 
  m_c[0]=v2(0); m_c[1]=v1(1); 
  m_h=v1.H()|v2.H(); m_s=v1.S()|v2.S();
}

template<class Scalar>
CAsT4<Scalar>::CAsT4(const CAsT4 &v,const SComplex &c)
{
  m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c;
  m_x[3]=v[3]*c; m_x[4]=v[4]*c; m_x[5]=v[5]*c;
  m_c[0]=v(0); m_c[1]=v(1); m_h=v.m_h; m_s=v.m_s;
}

template<class Scalar>
CAsT4<Scalar>::CAsT4(const CAsT4 &v,const Scalar &c)
{
  m_x[0]=v[0]*c; m_x[1]=v[1]*c; m_x[2]=v[2]*c;
  m_x[3]=v[3]*c; m_x[4]=v[4]*c; m_x[5]=v[5]*c;
  m_c[0]=v(0); m_c[1]=v(1); m_h=v.m_h; m_s=v.m_s;
}

template <class Scalar>
CAsT4<Scalar>& CAsT4<Scalar>::operator*=(const SComplex &c) 
{
  m_x[0]*=c; m_x[1]*=c; m_x[2]*=c; m_x[3]*=c; m_x[4]*=c; m_x[5]*=c;
  return *this;
}

template <class Scalar> CVec4<Scalar>
METOOLS::operator*(const CVec4<Scalar> &v,const CAsT4<Scalar> &t)
{ 
  CVec4<Scalar> j;
  j[0]=-t[0]/*-0,1*/*v[1]-t[3]/*-0,2*/*v[2]-t[5]/*-0,3*/*v[3];
  j[1]=-t[0]/*-0,1*/*v[0]-t[1]/*-1,2*/*v[2]-t[4]/*-1,3*/*v[3];
  j[2]=-t[3]/*-0,2*/*v[0]+t[1]/* 1,2*/*v[1]-t[2]/*-2,3*/*v[3];
  j[3]=-t[5]/*-0,3*/*v[0]+t[4]/* 1,3*/*v[1]+t[2]/* 2,3*/*v[2];
  return j;
}

template <class Scalar>
void CAsT4<Scalar>::Add(const CObject *c)
{
  const CAsT4 *t(static_cast<const CAsT4*>(c));
  m_x[0]+=t->m_x[0]; 
  m_x[1]+=t->m_x[1];
  m_x[2]+=t->m_x[2]; 
  m_x[3]+=t->m_x[3];
  m_x[4]+=t->m_x[4];
  m_x[5]+=t->m_x[5];
}

template <class Scalar>
void CAsT4<Scalar>::Divide(const double &d)
{
  m_x[0]/=Scalar(d);
  m_x[1]/=Scalar(d); 
  m_x[2]/=Scalar(d); 
  m_x[3]/=Scalar(d);
  m_x[4]/=Scalar(d);
  m_x[5]/=Scalar(d);
}

template <class Scalar>
void CAsT4<Scalar>::Multiply(const Complex &c)
{
  m_x[0]*=SComplex(c);
  m_x[1]*=SComplex(c); 
  m_x[2]*=SComplex(c); 
  m_x[3]*=SComplex(c);
  m_x[4]*=SComplex(c);
  m_x[5]*=SComplex(c);
}

template <class Scalar>
void CAsT4<Scalar>::Invert()
{
  m_x[0]=-m_x[0]; 
  m_x[1]=-m_x[1]; 
  m_x[2]=-m_x[2]; 
  m_x[3]=-m_x[3]; 
  m_x[4]=-m_x[4]; 
  m_x[5]=-m_x[5]; 
}

template <class Scalar>
bool CAsT4<Scalar>::IsZero() const
{
  return m_x[0]==Scalar(0.0) &&
    m_x[1]==Scalar(0.0) &&
    m_x[2]==Scalar(0.0) &&
    m_x[3]==Scalar(0.0) &&
    m_x[4]==Scalar(0.0) &&
    m_x[5]==Scalar(0.0);
}

template <class Scalar>
typename ATOOLS::AutoDelete_Vector<CAsT4<Scalar> > 
CAsT4<Scalar>::s_objects;

template <class Scalar>
CAsT4<Scalar> *CAsT4<Scalar>::New()
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new CAsT4();
  CAsT4 *v(s_objects.back());
  s_objects.pop_back();
  return v;
}

template <class Scalar>
CAsT4<Scalar> *CAsT4<Scalar>::New(const CAsT4 &s)
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new CAsT4(s);
  CAsT4 *v(s_objects.back());
  s_objects.pop_back();
  *v=s;
  return v;
}

template <class Scalar>
CObject *CAsT4<Scalar>::Copy() const
{
  return CAsT4::New(*this);
}

template <class Scalar>
void CAsT4<Scalar>::Delete()
{
#ifndef USING__Threading
  s_objects.push_back(this);
#else
  delete this;
#endif
}

namespace METOOLS {

  template class DCAsT4D;
  template CVec4<double> operator*
  (const CVec4<double> &v,const CAsT4<double> &t);
  template std::ostream &operator<<(std::ostream &ostr,const DCAsT4D &s);

  template class QCAsT4D;
  template CVec4<long double> operator*
  (const CVec4<long double> &v,const CAsT4<long double> &t);
  template std::ostream &operator<<(std::ostream &ostr,const QCAsT4D &s);

}
