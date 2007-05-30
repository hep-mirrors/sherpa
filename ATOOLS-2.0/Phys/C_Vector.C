#include "C_Vector.H"

using namespace ATOOLS;

double CVec4D::s_accu=1.0e-12;

std::ostream &ATOOLS::operator<<(std::ostream &s,const CVec4D &vec)
{
  return s<<'{'<<vec.H(0)<<","<<vec.H(1)<<";"<<vec(0)<<","<<vec(1)<<'|'
	  <<vec[0]<<','<<vec[1]<<','<<vec[2]<<','<<vec[3]<<'}';
}

bool CVec4D::IsZero() const
{
  for(short unsigned int i(0);i<4;++i) 
    if (std::abs(m_x[i])>Accu()) return false;
  return true;
}

bool CVec4D::Nan() const
{
  for(short unsigned int i(0);i<4;++i) {
    if (!(m_x[i].real()>=0.0) && !(m_x[i].real()<=0.0)) return true;
    if (!(m_x[i].imag()>=0.0) && !(m_x[i].imag()<=0.0)) return true;
  }
  return false;
}

void CVec4D::ResetAccu()                
{ 
  s_accu=Accu(); 
}

