#include "HELICITIES/Main/C_Vector.H"

using namespace ATOOLS;

template <class Scalar>
double CVec4<Scalar>::s_accu(1.0e-12);

template <class Scalar> std::ostream &
ATOOLS::operator<<(std::ostream &s,const CVec4<Scalar> &vec)
{
  return s<<'{'<<vec.H(0)<<","<<vec.H(1)<<";"<<vec(0)<<","<<vec(1)<<'|'
	  <<vec[0]<<','<<vec[1]<<','<<vec[2]<<','<<vec[3]<<'}';
}

template <class Scalar>
bool CVec4<Scalar>::IsZero() const
{
  for(short unsigned int i(0);i<4;++i) 
    if (!ATOOLS::IsZero(m_x[i])) return false;
  return true;
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

namespace ATOOLS {

  template class DCVec4D;
  template std::ostream &operator<<(std::ostream &ostr,const DCVec4D &s);

  template class QCVec4D;
  template std::ostream &operator<<(std::ostream &ostr,const QCVec4D &s);

}
