#include "SINIC++/Cross_Sections/Sigma_Total.H"
#include "SINIC++/Tools/Special_Functions.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SINIC;
using namespace ATOOLS;

double Sigma_Tot::GetValue(const double & B) { 
  return 2.*p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B)/2.)); 
}

double Sigma_Tot::GetCombinedValue(const double & B) { 
  double value(0.); //, pref, eik;
  for (std::list<Omega_ik *>::iterator eikonal=p_eikonals->begin();
       eikonal!=p_eikonals->end(); eikonal++) {
    value += 2.*(*eikonal)->Prefactor()*(1.-exp(-(**eikonal)(B)/2.)); 
  }
  return value;
}

void Sigma_Tot::TestTotalCrossSection(){
  const double EulerGamma= 0.577215664901532860606512090082 ;
  double m_a,m_c,m_alpha,m_res,m_ei;
  double m_Delta,m_prefactor,m_Lambda2,m_beta0,m_kappa;
  ExpInt m_expint;
  m_Delta     = (*p_eikonals).front()->Delta();
  m_prefactor = (*p_eikonals).front()->Prefactor();
  m_kappa     = (*p_eikonals).front()->Kappa_i();
  m_Lambda2   = (*p_eikonals).front()->Lambda2();
  m_beta0     = (*p_eikonals).front()->FF1()->Beta0();
  m_a         = m_Lambda2/(8.*(1.+m_kappa));
  m_c         =
    ATOOLS::sqr(m_beta0)*m_Lambda2*(1.+m_kappa)*exp(2.*m_Delta*m_Y)/(8.*M_PI);
  m_alpha     = 4.*M_PI*m_prefactor;
  m_ei        = m_expint.GetExpInt(-m_c/2.);
  m_res       = m_alpha*(EulerGamma-m_ei+log(m_c/2.))/(2.*m_a);
  msg_Out() << "In " << METHOD << " sigma_tot = "<< m_res <<" 1/GeV^2 = "
	    <<m_res*rpa->Picobarn()/1.e9<<" mb ."<<std::endl;
  //	exit(1);
}


double Elastic_Slope::GetValue(const double & B) { 
  return B*B*p_eikonal->Prefactor()*(1.-exp(-(*p_eikonal)(B)/2.))/m_stot; 
}

double Elastic_Slope::GetCombinedValue(const double & B) { 
  double value(0.); //, pref, eik;
  for (std::list<Omega_ik *>::iterator eikonal=p_eikonals->begin();
       eikonal!=p_eikonals->end(); eikonal++) {
    value += (*eikonal)->Prefactor()*(1.-exp(-(**eikonal)(B)/2.)); 
  }
  return B*B*value/m_stot;
}

