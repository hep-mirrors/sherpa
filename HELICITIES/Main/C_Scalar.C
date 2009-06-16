#include "HELICITIES/Main/C_Scalar.H"

using namespace ATOOLS;

template <class Scalar>
double CScalar<Scalar>::s_accu(1.0e-12);

template <class Scalar> std::ostream &
ATOOLS::operator<<(std::ostream &str,const CScalar<Scalar> &s)
{
  return str<<'{'<<s.H(0)<<","<<s.H(1)<<'|'<<s[0]<<'}';
}

template <class Scalar>
void CScalar<Scalar>::ResetAccu()                
{ 
  s_accu=Accu(); 
}

namespace ATOOLS {

  template class DCScalar;
  template std::ostream &operator<<(std::ostream &ostr,const DCScalar &s);

  template class QCScalar;
  template std::ostream &operator<<(std::ostream &ostr,const QCScalar &s);

}
