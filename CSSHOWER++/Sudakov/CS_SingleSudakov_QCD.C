#include "CSSHOWER++/Sudakov/CS_SingleSudakov_QCD.H"

#include "ATOOLS/Math/Gauss_Integrator.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

//###########################################################################
//###########################################################################
//the final state emitter -- final state spectator splittings 
//###########################################################################
//###########################################################################

//###########################################################################
// q->qg
//###########################################################################
DeltaFF_q_qg::DeltaFF_q_qg(const ATOOLS::Flavour & flav) : m_gauss(this)
{
  m_type = cstp::FF;
  m_mi2  = sqr(flav.Mass());
  m_fl   = flav; 
}

double DeltaFF_q_qg::operator()(double kt2min,double kt2max) {

  return 1.;
}

double DeltaFF_q_qg::operator()(double kt2min) {
  return operator()(kt2min,m_kt2max);
}

double DeltaFF_q_qg::Delta(double kt2min,double kt2max) {
  SetKT2max(kt2max);
  double integral = m_gauss.Integrate(kt2max,kt2min,1.e-4,1);
  return exp(-integral);  
}

//###########################################################################
// g->gg
//###########################################################################
DeltaFF_g_gg::DeltaFF_g_gg() : m_gauss(this)
{
  m_type = cstp::FF;
  m_fl   = Flavour(kf_gluon); 
}

double DeltaFF_g_gg::operator()(double kt2min,double kt2max) {
  return 1.;
}

double DeltaFF_g_gg::operator()(double kt2min) {
  return operator()(kt2min,m_kt2max);
}
  
double DeltaFF_g_gg::Delta(double kt2min,double kt2max) {
  SetKT2max(kt2max);
  double integral = m_gauss.Integrate(kt2max,kt2min,1.e-4,1);
  return exp(-integral);  
}


//###########################################################################
// g->qq
//###########################################################################
DeltaFF_g_qq::DeltaFF_g_qq(const ATOOLS::Flavour & flav) : m_gauss(this) 
{
  m_type = cstp::FF;
  m_mj2  = sqr(flav.Mass());
  m_fl   = Flavour(kf_gluon); 
}

double DeltaFF_g_qq::operator()(double kt2min,double kt2max) {
  
  //std::cout<<" called operator ... "<<std::endl; 
  double kt2max_local=kt2max;
  if (kt2max>m_Q2) {
    kt2max_local=m_Q2/2.;
  }
  double zm = 0.5*(1.-sqrt(1.-kt2min/kt2max_local));
  double zp = 0.5*(1.+sqrt(1.-kt2min/kt2max_local));
  double G_zintegrated = s_TR/(2.*M_PI)*(-zm+zp)/(3.*m_Q2*(zm-zp))*
    ((zm-zp)*(6.*kt2max_local+m_Q2*(3.+2.*zm*zm+(zm+zp)*(2.*zp-3.)))
     +3.*kt2max_local*log(((-1.+zm)*zp)/(zm*(-1.+zp))));
  
  if (G_zintegrated<0. || (G_zintegrated<0. && G_zintegrated>0.)) {
    //std::cout<<" z's "<<zm<<" "<<zp<<std::endl;  
    //std::cout<<" Gamma_g_qq is : "<<G_zintegrated<<std::endl; 
    G_zintegrated=0.;
  }
  return G_zintegrated;
}

double DeltaFF_g_qq::operator()(double kt2min) {
  return operator()(kt2min,m_kt2max);
}
  
double DeltaFF_g_qq::Delta(double kt2min,double kt2max) {
  
  SetKT2max(kt2max);
  double integral = m_gauss.Integrate(kt2min,kt2max,1.e-3,1);
  //std::cout<<" integral yields : "<<integral<<std::endl; 
  return exp(-integral);  
}
