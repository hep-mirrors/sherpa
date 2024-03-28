#include "YFS/Main/Weight.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Particle.H"








Weight::Weight(){

}


Weight::~Weight(){

}


voif Weight::TotalWeight(double v, double g, double gp, double eps){
  double expF = exp(m_g*(1.0/4.) + m_alp/M_PI*(pow(M_PI,2.)/3.-0.5));
  if(m_v > m_eps){
    m_weight=m_gp*pow(m_v,m_gp-1.)*diljac0*pow(m_eps,m_g-m_gp);
  }
  else{
    m_massW = 1.0;
    jacW = 1.0;
    m_weight =  m_g*pow(m_v,m_g-1.);
  }
  m_weight*=expF*m_cut*m_massW*jacW;
  DEBUG_FUNC("Born Weight = "<<m_born<<std::endl<<
       "v = "<<m_v<<std::endl<<
       "vmin = "<<m_vmin<<std::endl<<
       "vmax = "<<m_vmax<<std::endl<<
       "eps  = "<<m_vmin<<std::endl<<
       "Kmin = "<<m_Kmin<<std::endl<<
       "Lambda = "<<lam<<std::endl << "J0 = "<<diljac0<<std::endl<<
       "J = "<<diljac<<std::endl << "W_J = "<<jacW<<std::endl<<"cut = "
       <<m_cut<<std::endl
       << "W_mass = "<<m_massW<<std::endl);
  for (int i = 0; i < m_photons.size(); ++i) msg_Debugging()<<"k["<< i<<"] = "<<m_photons[i]<<std::endl;
  if(IsBad(m_weight)){
    msg_Error()<<"YFS Weight is: "<<m_weight<<std::endl<<
    "Born Weight = "<<m_born<<std::endl<<
       "v = "<<m_v<<std::endl<<
       "vmin = "<<m_vmin<<std::endl<<
       "vmax = "<<m_vmax<<std::endl<<
       "eps  = "<<m_vmin<<std::endl<<
       "Kmin = "<<m_Kmin<<std::endl<<
       "Lambda = "<<lam<<std::endl << "J0 = "<<diljac0<<std::endl<<
       "J = "<<diljac<<std::endl << "W_J = "<<jacW<<std::endl<<"cut = "
       <<m_cut<<std::endl
       << "W_mass = "<<m_massW<<std::endl;
  }
}