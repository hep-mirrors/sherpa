#include "METOOLS/Currents/C_Vector.H"

#include "ATOOLS/Org/CXXFLAGS.H"

using namespace METOOLS;

template <class Scalar>
double CVec4<Scalar>::s_accu(1.0e-12);

template <class Scalar> std::ostream &
METOOLS::operator<<(std::ostream &s,const CVec4<Scalar> &vec)
{
  return s<<'{'<<vec.H()<<","<<vec.S()<<";"<<vec(0)<<","<<vec(1)<<'|'
	  <<vec[0]<<','<<vec[1]<<','<<vec[2]<<','<<vec[3]<<'}';
}

template <class Scalar>
bool CVec4<Scalar>::Nan() const
{
  for(short unsigned int i(0);i<4;++i) {
    if (ATOOLS::IsNan(m_x[i])) return true;
  }
  return false;
}

template <class Scalar>
void CVec4<Scalar>::ResetAccu()                
{ 
  s_accu=Accu(); 
}

template <class Scalar>
void CVec4<Scalar>::Add(const CObject *c)
{
  const CVec4 *v(static_cast<const CVec4*>(c));
  m_x[0]+=v->m_x[0]; 
  m_x[1]+=v->m_x[1];
  m_x[2]+=v->m_x[2]; 
  m_x[3]+=v->m_x[3];
}

template <class Scalar>
void CVec4<Scalar>::Divide(const double &d)
{
  m_x[0]/=Scalar(d);
  m_x[1]/=Scalar(d); 
  m_x[2]/=Scalar(d); 
  m_x[3]/=Scalar(d);
}

template <class Scalar>
void CVec4<Scalar>::Multiply(const Complex &c)
{
  m_x[0]*=SComplex(c);
  m_x[1]*=SComplex(c); 
  m_x[2]*=SComplex(c); 
  m_x[3]*=SComplex(c);
}

template <class Scalar>
void CVec4<Scalar>::Invert()
{
  m_x[0]=-m_x[0]; 
  m_x[1]=-m_x[1]; 
  m_x[2]=-m_x[2]; 
  m_x[3]=-m_x[3]; 
}

template <class Scalar>
bool CVec4<Scalar>::IsZero() const
{
  return m_x[0]==Scalar(0.0) &&
    m_x[1]==Scalar(0.0) &&
    m_x[2]==Scalar(0.0) &&
    m_x[3]==Scalar(0.0);
}

template <class Scalar>
typename ATOOLS::AutoDelete_Vector<CVec4<Scalar> > 
CVec4<Scalar>::s_objects;

template <class Scalar>
CVec4<Scalar> *CVec4<Scalar>::New()
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new CVec4();
  CVec4 *v(s_objects.back());
  s_objects.pop_back();
  return v;
}

template <class Scalar>
CVec4<Scalar> *CVec4<Scalar>::New(const CVec4 &s)
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new CVec4(s);
  CVec4 *v(s_objects.back());
  s_objects.pop_back();
  *v=s;
  return v;
}

template <class Scalar>
CVec4<Scalar> *CVec4<Scalar>::New
(const Scalar &x0, const Scalar &x1, 
 const Scalar &x2, const Scalar &x3,
 const int c1,const int c2,
 const size_t &h,const size_t &s)
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new CVec4(x0,x1,x2,x3,c1,c2,h,s);
  CVec4 *v(s_objects.back());
  s_objects.pop_back();
  v->m_x[0]=x0;
  v->m_x[1]=x1;
  v->m_x[2]=x2;
  v->m_x[3]=x3; 
  v->m_c[0]=c1;
  v->m_c[1]=c2;
  v->m_h=h;
  v->m_s=s;
  return v;
}

template <class Scalar>
CObject *CVec4<Scalar>::Copy() const
{
  return CVec4::New(*this);
}

template <class Scalar>
void CVec4<Scalar>::Delete()
{
#ifndef USING__Threading
  s_objects.push_back(this);
#else
  delete this;
#endif
}

namespace METOOLS {

  template class DCVec4D;
  template std::ostream &operator<<(std::ostream &ostr,const DCVec4D &s);

  template class QCVec4D;
  template std::ostream &operator<<(std::ostream &ostr,const QCVec4D &s);

}
