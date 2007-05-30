#include "C_Tensor.H"

using namespace ATOOLS;

double CAsT4D::s_accu=1.0e-12;

std::ostream &ATOOLS::operator<<(std::ostream &s,const CAsT4D &ten)
{
  return s<<'{'<<ten.H(0)<<","<<ten.H(1)<<";"<<ten(0)<<","<<ten(1)<<'|'
	  <<ten[0]<<','<<ten[1]<<','<<ten[2]<<','
	  <<ten[3]<<','<<ten[4]<<','<<ten[5]<<'}';
}

bool CAsT4D::IsZero() const
{
  for(short unsigned int i(0);i<6;++i) 
    if (std::abs(m_x[i])>Accu()) return false;
  return true;
}

bool CAsT4D::Nan() const
{
  for(short unsigned int i(0);i<6;++i) {
    if (!(m_x[i].real()>=0.0) && !(m_x[i].real()<=0.0)) return true;
    if (!(m_x[i].imag()>=0.0) && !(m_x[i].imag()<=0.0)) return true;
  }
  return false;
}

void CAsT4D::ResetAccu()                
{ 
  s_accu=Accu(); 
}

