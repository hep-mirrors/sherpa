#include "HELICITIES/Main/C_Tensor.H"

using namespace ATOOLS;

template<class Scalar>
double CAsT4<Scalar>::s_accu(1.0e-12);

template<class Scalar> std::ostream &
ATOOLS::operator<<(std::ostream &s,const CAsT4<Scalar> &ten)
{
  return s<<'{'<<ten.H(0)<<","<<ten.H(1)<<";"<<ten(0)<<","<<ten(1)<<'|'
	  <<ten[0]<<','<<ten[1]<<','<<ten[2]<<','
	  <<ten[3]<<','<<ten[4]<<','<<ten[5]<<'}';
}

template<class Scalar>
bool CAsT4<Scalar>::IsZero() const
{
  for(short unsigned int i(0);i<6;++i) 
    if (std::abs(m_x[i])>Accu()) return false;
  return true;
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

namespace ATOOLS {

  template class DCAsT4D;
  template std::ostream &operator<<(std::ostream &ostr,const DCAsT4D &s);

  template class QCAsT4D;
  template std::ostream &operator<<(std::ostream &ostr,const QCAsT4D &s);

}
