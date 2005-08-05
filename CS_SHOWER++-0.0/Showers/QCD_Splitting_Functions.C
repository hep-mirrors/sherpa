#include "QCD_Splitting_Functions.H"
#include "Random.H"


using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

//###########################################################################
// q->qg
//###########################################################################
q_qg_FF::q_qg_FF(const ATOOLS::Flavour & flav) : m_Jmax(0.5)
{
  m_type     = cstp::FF;
  m_flavs[0] = flav;
  m_flavs[1] = flav;
  m_flavs[2] = Flavour(kf::gluon);
}

double q_qg_FF::operator() (const double z,const double y) {
  return s_CF* ( 2./(1.-z+z*y) - (1.+z) ) * J(y);
}

double q_qg_FF::OverIntegrated(const double zmin,const double zmax,const double scale) {
  m_lastint = m_zmin = m_zmax = 0.; 
  //if (scale<m_flavs[0].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_lastint = 2.*s_CF*log((1.-zmin)/(1.-zmax)) * m_Jmax;
  return m_lastint;
}

double q_qg_FF::Overestimated(const double z,const double y) {
  return s_CF* ( 2./(1.-z) ) * m_Jmax;
}

double q_qg_FF::RejectionWeight(const double z,const double y) {
  return ( (1.-z)/(1.-z+z*y) - (1.-z*z)/2. ) * J(y)/m_Jmax;
}

double q_qg_FF::Z() {
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran.Get());
}

double q_qg_FF::J(const double y) { return m_Jmax*(1.-y); }


//###########################################################################
// q->gq
//###########################################################################
q_gq_FF::q_gq_FF(const ATOOLS::Flavour & flav) : m_Jmax(0.5)
{
  m_type     = cstp::FF;
  m_flavs[0] = flav;
  m_flavs[1] = Flavour(kf::gluon);
  m_flavs[2] = flav;
}

double q_gq_FF::operator() (const double z,const double y) {
  return s_CF* ( 2./(z+(1.-z)*y) - z ) * J(y);
}

double q_gq_FF::OverIntegrated(const double zmin,const double zmax,const double scale) {
  m_lastint = m_zmin = m_zmax = 0.; 
  //if (scale<m_flavs[0].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_lastint = 2.*s_CF*log(zmax/zmin) * m_Jmax;
  return m_lastint;
}

double q_gq_FF::Overestimated(const double z,const double y) {
  return s_CF* ( 2./z ) * m_Jmax;
}

double q_gq_FF::RejectionWeight(const double z,const double y) {
  return ( z/(z+(1.-z)*y) - z*z/2. ) * J(y)/m_Jmax;
}

double q_gq_FF::Z() {
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double q_gq_FF::J(const double y) { return m_Jmax*(1.-y); }


//###########################################################################
// g->gg
//###########################################################################
g_gg_FF::g_gg_FF() : m_Jmax(0.5)
{
  m_type     = cstp::FF;
  m_flavs[0] = Flavour(kf::gluon);
  m_flavs[1] = Flavour(kf::gluon);
  m_flavs[2] = Flavour(kf::gluon);
}

double g_gg_FF::operator() (const double z,const double y) {
  return s_CA * ( 1./(1.-z+z*y) + 1./(z+y-z*y) -2. + z*(1.-z) ) * J(y);
}

double g_gg_FF::OverIntegrated(const double zmin,const double zmax,const double scale) {
  m_zmin = zmin; m_zmax = zmax;
  m_lastint = s_CA * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
  return m_lastint;
}

double g_gg_FF::Overestimated(const double z,const double y) {
  return s_CA * ( 1./(1.-z) + 1./z ) * m_Jmax;
}

double g_gg_FF::RejectionWeight(const double z,const double y) {
  return z*(1.-z) * ( 1./(1.-z+z*y) + 1./(z+y-z*y) -2. + z*(1.-z) ) * J(y)/m_Jmax;
}

double g_gg_FF::Z() {
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double g_gg_FF::J(const double y) { return m_Jmax*(1.-y); }


//###########################################################################
// g->qq
//###########################################################################
g_qq_FF::g_qq_FF(const ATOOLS::Flavour & flav) : m_Jmax(0.5) 
{
  m_type     = cstp::FF;
  m_flavs[0] = Flavour(kf::gluon);
  m_flavs[1] = flav;
  m_flavs[2] = flav.Bar();
}

double g_qq_FF::operator() (const double z,const double y) {
  return s_TR * (1.-2.*z*(1.-z)) * J(y);
}

double g_qq_FF::OverIntegrated(const double zmin,const double zmax,const double scale) {
  m_lastint = m_zmin = m_zmax = 0.; 
  if (scale<m_flavs[1].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_lastint = s_TR*(m_zmax-m_zmin) * m_Jmax;                                             
  return m_lastint;
}

double g_qq_FF::Overestimated(const double z,const double y) {
  return s_TR * m_Jmax;
}

double g_qq_FF::RejectionWeight(const double z,const double y) {
  return (1.-2.*z*(1.-z)) * J(y)/m_Jmax;
}

double g_qq_FF::Z() {
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double g_qq_FF::J(const double y) { return m_Jmax*(1.-y); }
