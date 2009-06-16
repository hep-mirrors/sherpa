#include "AMISIC++/Model/Reggeon_Trajectory.H"

#include "ATOOLS/Math/MathTools.H"

using namespace AMISIC;

Reggeon_Trajectory::Reggeon_Trajectory(const double alpha_0,
				       const double alpha_prime):
  m_alpha_0(alpha_0), 
  m_alpha_prime(alpha_prime), 
  m_s_0(1.0) {}
    
Reggeon_Trajectory::~Reggeon_Trajectory() 
{
}
    
void Reggeon_Trajectory::Fit(const double t,const double dsigma)
{ 
  m_s_0=m_s/pow(dsigma,1.0/(2.0*(m_alpha_0-1.0-m_alpha_prime*t))); 
}
    
double Reggeon_Trajectory::DSigma(const double t) const
{ 
  return pow(m_s/m_s_0,2.0*(m_alpha_0-1.0-m_alpha_prime*t)); 
}
    
double Reggeon_Trajectory::GetT(const double tmin,const double tmax,
				const double ran) const
{
  double C=m_s/m_s_0, D=m_alpha_0-1.0;
  return (D-log(pow(C,2.0*(D-m_alpha_prime*tmin))*(1.0-ran)
		+pow(C,2.0*(D-m_alpha_prime*tmax))*ran)/(2.0*log(C)))/
    m_alpha_prime;
}
    
