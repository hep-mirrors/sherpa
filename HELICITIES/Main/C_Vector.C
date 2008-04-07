#include "C_Vector.H"

using namespace HELICITIES;

double CVec4D::s_accu=1.0e-12;

std::ostream &HELICITIES::operator<<(std::ostream &s,const CVec4D &vec)
{
  if(vec(0)==-1 && vec(1)==-1) {
    return s<<'('<<vec[0]<<','<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
  }
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

CVec4D HELICITIES::cross(const CVec4D& v1, const CVec4D& v2, const CVec4D& v3)
{
  // \epsilon^{\mu\nu\rho\sigma}  v1_\nu  v2_\rho  v3_\sigma
  Complex x0 = -v1[1]*v2[2]*v3[3] -v1[2]*v2[3]*v3[1] -v1[3]*v2[1]*v3[2] +v1[1]*v2[3]*v3[2] +v1[3]*v2[2]*v3[1] +v1[2]*v2[1]*v3[3];
  Complex x1 = -v1[0]*v2[2]*v3[3] -v1[2]*v2[3]*v3[0] -v1[3]*v2[0]*v3[2] +v1[0]*v2[3]*v3[2] +v1[3]*v2[2]*v3[0] +v1[2]*v2[0]*v3[3];
  Complex x2 =  v1[0]*v2[1]*v3[3] +v1[1]*v2[3]*v3[0] +v1[3]*v2[0]*v3[1] -v1[0]*v2[3]*v3[1] -v1[3]*v2[1]*v3[0] -v1[1]*v2[0]*v3[3];
  Complex x3 = -v1[0]*v2[1]*v3[2] -v1[1]*v2[2]*v3[0] -v1[2]*v2[0]*v3[1] +v1[0]*v2[2]*v3[1] +v1[2]*v2[1]*v3[0] +v1[1]*v2[0]*v3[2];
  return CVec4D(x0,x1,x2,x3);
}

