#include "METOOLS/Currents/C_Scalar.H"

using namespace METOOLS;

template <class Scalar>
double CScalar<Scalar>::s_accu(1.0e-12);

template <class Scalar> std::ostream &
METOOLS::operator<<(std::ostream &str,const CScalar<Scalar> &s)
{
  return str<<'{'<<s.H()<<','<<s.S()<<';'<<s(0)<<','<<s(1)<<'|'<<s[0]<<'}';
}

template <class Scalar>
void CScalar<Scalar>::ResetAccu()                
{ 
  s_accu=Accu(); 
}

template <class Scalar>
void CScalar<Scalar>::Add(const CObject *c)
{
  const CScalar *s(static_cast<const CScalar*>(c));
  m_x+=s->m_x;
}

template <class Scalar>
void CScalar<Scalar>::Divide(const double &d)
{
  m_x/=Scalar(d);
}

template <class Scalar>
void CScalar<Scalar>::Invert()
{
  m_x=-m_x;
}

template <class Scalar>
typename ATOOLS::AutoDelete_Vector<CScalar<Scalar> > 
CScalar<Scalar>::s_objects;

template <class Scalar>
bool CScalar<Scalar>::IsZero() const
{
  return m_x==Scalar(0.0);
}

template <class Scalar>
CScalar<Scalar> *CScalar<Scalar>::New()
{
  s_objects.MtxLock();
  if (s_objects.empty()) {
    s_objects.MtxUnLock();
    return new CScalar();
  }
  CScalar *v(s_objects.back());
  s_objects.pop_back();
  s_objects.MtxUnLock();
  return v;
}

template <class Scalar>
CScalar<Scalar> *CScalar<Scalar>::New(const CScalar &s)
{
  s_objects.MtxLock();
  if (s_objects.empty()) {
    s_objects.MtxUnLock();
    return new CScalar(s);
  }
  CScalar *v(s_objects.back());
  s_objects.pop_back();
  s_objects.MtxUnLock();
  *v=s;
  return v;
}

template <class Scalar>
CObject *CScalar<Scalar>::Copy() const
{
  return CScalar::New(*this);
}

template <class Scalar>
void CScalar<Scalar>::Delete()
{
  s_objects.MtxLock();
  s_objects.push_back(this);
  s_objects.MtxUnLock();
}

namespace METOOLS {

  template class DCScalar;
  template std::ostream &operator<<(std::ostream &ostr,const DCScalar &s);

  template class QCScalar;
  template std::ostream &operator<<(std::ostream &ostr,const QCScalar &s);

}
