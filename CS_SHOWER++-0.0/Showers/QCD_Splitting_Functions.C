#include "QCD_Splitting_Functions.H"
#include "Random.H"


using namespace CS_SHOWER;
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
q_qg_FF::q_qg_FF(const ATOOLS::Flavour & flav) : m_Jmax(0.5)
{
  m_type     = cstp::FF;
  m_flavs[0] = flav;
  m_flavs[1] = flav;
  m_flavs[2] = Flavour(kf::gluon);
}

double q_qg_FF::operator() (const double z,const double y,
			    const double eta,const double scale) {
  return s_CF* ( 2./(1.-z+z*y) - (1.+z) ) * J(y);
}

double q_qg_FF::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  //if (scale<m_flavs[0].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_lastint = 2.*s_CF*log((1.-zmin)/(1.-zmax)) * m_Jmax;
  return m_lastint;
}

double q_qg_FF::Overestimated(const double z,const double y) {
  return s_CF* ( 2./(1.-z) ) * m_Jmax;
}

double q_qg_FF::RejectionWeight(const double z,const double y,
				const double eta,const double scale) {
  return ( (1.-z)/(1.-z+z*y) - (1.-z*z)/2. ) * J(y)/m_Jmax;
}

double q_qg_FF::Z() {
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran.Get());
}

double q_qg_FF::J(const double y,const double eta,const double scale) { return m_Jmax*(1.-y); }


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

double q_gq_FF::operator() (const double z,const double y,const double eta,const double scale) {
  return s_CF* ( 2./z - 2 + z ) * J(y);
}

double q_gq_FF::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  //if (scale<m_flavs[0].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_lastint = 2.*s_CF*log(zmax/zmin) * m_Jmax;
  return m_lastint;
}

double q_gq_FF::Overestimated(const double z,const double y) {
  return s_CF* ( 2./z ) * m_Jmax;
}

double q_gq_FF::RejectionWeight(const double z,const double y,
				const double eta,const double scale) {
  return ( 1. - z + z*z/2. ) * J(y)/m_Jmax;
}

double q_gq_FF::Z() {
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double q_gq_FF::J(const double y,const double eta,const double scale) { return m_Jmax*(1.-y); }


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

double g_gg_FF::operator() (const double z,const double y,const double eta,const double scale) {
  return s_CA * ( 1./(1.-z+z*y) + 1./(z+y-z*y) -2. + z*(1.-z) ) * J(y);
}

double g_gg_FF::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_zmin = zmin; m_zmax = zmax;
  m_lastint = s_CA * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
  return m_lastint;
}

double g_gg_FF::Overestimated(const double z,const double y) {
  return s_CA * ( 1./(1.-z) + 1./z ) * m_Jmax;
}

double g_gg_FF::RejectionWeight(const double z,const double y,
				const double eta,const double scale) {
  return z*(1.-z) * ( 1./(1.-z+z*y) + 1./(z+y-z*y) -2. + z*(1.-z) ) * J(y)/m_Jmax;
}

double g_gg_FF::Z() {
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double g_gg_FF::J(const double y,const double eta,const double scale) { return m_Jmax*(1.-y); }


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

double g_qq_FF::operator() (const double z,const double y,
			    const double eta,const double scale) {
  return s_TR * (1.-2.*z*(1.-z)) * J(y);
}

double g_qq_FF::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  if (scale<m_flavs[1].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_lastint = s_TR*(m_zmax-m_zmin) * m_Jmax;                                             
  return m_lastint;
}

double g_qq_FF::Overestimated(const double z,const double y) {
  return s_TR * m_Jmax;
}

double g_qq_FF::RejectionWeight(const double z,const double y,
				const double eta,const double scale) {
  return (1.-2.*z*(1.-z)) * J(y)/m_Jmax;
}

double g_qq_FF::Z() {
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double g_qq_FF::J(const double y,const double eta,const double scale) { return m_Jmax*(1.-y); }

//###########################################################################
//###########################################################################
//the final state emitter -- initial state spectator splittings 
//###########################################################################
//###########################################################################

//###########################################################################
// q->qg
//###########################################################################
q_qg_FI::q_qg_FI(const ATOOLS::Flavour & flav,PDF::PDF_Base * pdf) : 
  m_Jmax(1.), p_pdf(pdf)
{
  m_type     = cstp::FI;
  m_flavs[0] = flav;
  m_flavs[1] = flav;
  m_flavs[2] = Flavour(kf::gluon);
}

double q_qg_FI::operator() (const double z,const double y,
			    const double eta, const double scale) {
  return s_CF* ( 2./(1.-z+y) - (1.+z) ) * J(y,eta,scale);
}

double q_qg_FI::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1./(2.*xbj); 
  m_lastint = 2.*s_CF*log((1.-zmin)/(1.-zmax)) * m_Jmax;
  return m_lastint;
}

double q_qg_FI::Overestimated(const double z,const double y) {
  return s_CF* ( 2./(1.-z) ) * m_Jmax;
}

double q_qg_FI::RejectionWeight(const double z,const double y,
				const double eta, const double scale) {
  return ( (1.-z)/(1.-z+y) - (1.-z*z)/2. ) * J(y,eta,scale)/m_Jmax;
}

double q_qg_FI::Z() {
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran.Get());
}

double q_qg_FI::J(const double y,const double eta,const double scale) { 
  
  p_pdf->Calculate(eta/(1.-y),scale);
  double fresh = p_pdf->GetXPDF(m_flavs[0]);
  
  p_pdf->Calculate(eta,scale);
  double old = p_pdf->GetXPDF(m_flavs[0]);

  return 0.5/(1.-y) * fresh/old;
}

//###########################################################################
// q->gq
//###########################################################################
q_gq_FI::q_gq_FI(const ATOOLS::Flavour & flav,PDF::PDF_Base * pdf) : 
  m_Jmax(0.5), p_pdf(pdf)
{
  m_type     = cstp::FI;
  m_flavs[0] = flav;
  m_flavs[1] = Flavour(kf::gluon);
  m_flavs[2] = flav;
}

double q_gq_FI::operator() (const double z,const double y,
			    const double eta, const double scale) {
  return s_CF* ( 2./z - 2 + z ) * J(y,eta,scale);
}

double q_gq_FI::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  //if (scale<m_flavs[0].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1./(2.*xbj); 
  m_lastint = 2.*s_CF*log(zmax/zmin) * m_Jmax;
  return m_lastint;
}

double q_gq_FI::Overestimated(const double z,const double y) {
  return s_CF* ( 2./z ) * m_Jmax;
}

double q_gq_FI::RejectionWeight(const double z,const double y,
				const double eta, const double scale) {
  return ( 1. - z + z*z/2. ) * J(y,eta,scale)/m_Jmax;
}

double q_gq_FI::Z() {
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double q_gq_FI::J(const double y,const double eta,const double scale) {   
  p_pdf->Calculate(eta/(1.-y),scale);
  double fresh = p_pdf->GetXPDF(m_flavs[0]);
  
  p_pdf->Calculate(eta,scale);
  double old = p_pdf->GetXPDF(m_flavs[0]);

  return 0.5/(1.-y) * fresh/old;
}

//###########################################################################
// g->gg
//###########################################################################
g_gg_FI::g_gg_FI(PDF::PDF_Base * pdf) : m_Jmax(0.5), p_pdf(pdf)
{
  m_type     = cstp::FI;
  m_flavs[0] = Flavour(kf::gluon);
  m_flavs[1] = Flavour(kf::gluon);
  m_flavs[2] = Flavour(kf::gluon);
}

double g_gg_FI::operator() (const double z,const double y,
			    const double eta, const double scale) {
  return s_CA * ( 1./(1.-z+y) + 1./(z+y) -2. + z*(1.-z) ) * J(y,eta,scale);
}

double g_gg_FI::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1./(2.*xbj); 
  m_lastint = s_CA * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
  return m_lastint;
}

double g_gg_FI::Overestimated(const double z,const double y) {
  return s_CA * ( 1./(1.-z) + 1./z ) * m_Jmax;
}

double g_gg_FI::RejectionWeight(const double z,const double y,
				const double eta, const double scale) {
  return z*(1.-z) * ( 1./(1.-z+y) + 1./(z+y) -2. + z*(1.-z) ) * J(y,eta,scale)/m_Jmax;
}

double g_gg_FI::Z() {
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double g_gg_FI::J(const double y,const double eta,const double scale) { 
  
  p_pdf->Calculate(eta/(1.-y),scale);
  double fresh = p_pdf->GetXPDF(m_flavs[0]);
  
  p_pdf->Calculate(eta,scale);
  double old = p_pdf->GetXPDF(m_flavs[0]);

  return 0.5/(1.-y) * fresh/old;
}


//###########################################################################
// g->qq
//###########################################################################
g_qq_FI::g_qq_FI(const ATOOLS::Flavour & flav,PDF::PDF_Base * pdf) : 
  m_Jmax(0.5), p_pdf(pdf)
{
  m_type     = cstp::FI;
  m_flavs[0] = Flavour(kf::gluon);
  m_flavs[1] = flav;
  m_flavs[2] = flav.Bar();
}

double g_qq_FI::operator() (const double z,const double y,
			    const double eta, const double scale) {
  return s_TR * (1.-2.*z*(1.-z)) * J(y,eta,scale);
}

double g_qq_FI::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  if (scale<m_flavs[1].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1./(2.*xbj); 
  m_lastint = s_TR*(m_zmax-m_zmin) * m_Jmax;                                             
  return m_lastint;
}

double g_qq_FI::Overestimated(const double z,const double y) {
  return s_TR * m_Jmax;
}

double g_qq_FI::RejectionWeight(const double z,const double y,
				const double eta, const double scale) {
  return (1.-2.*z*(1.-z)) * J(y,eta,scale)/m_Jmax;
}

double g_qq_FI::Z() {
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double g_qq_FI::J(const double y,const double eta,const double scale) { 
  
  p_pdf->Calculate(eta/(1.-y),scale);
  double fresh = p_pdf->GetXPDF(m_flavs[0]);
  
  p_pdf->Calculate(eta,scale);
  double old = p_pdf->GetXPDF(m_flavs[0]);

  return 0.5/(1.-y) * fresh/old;
}

//###########################################################################
//###########################################################################
//the initial state emitter -- final state spectator splittings 
//###########################################################################
//###########################################################################

//###########################################################################
// q->qg
//###########################################################################
q_qg_IF::q_qg_IF(const ATOOLS::Flavour & flav,PDF::PDF_Base * pdf) : 
  m_Jmax(1.), p_pdf(pdf)
{
  m_type     = cstp::IF;
  m_flavs[0] = flav;
  m_flavs[1] = flav;
  m_flavs[2] = Flavour(kf::gluon);
}

double q_qg_IF::operator() (const double z,const double y,
			    const double eta, const double scale) {
  return s_CF* ( 2./(1.-z+y) - (1.+z) ) * J(z,eta,scale);
}

double q_qg_IF::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1./xbj; 
  m_lastint = 2.*s_CF*log((1.-zmin)/(1.-zmax)) * m_Jmax;
  return m_lastint;
}

double q_qg_IF::Overestimated(const double z,const double y) {
  return s_CF* ( 2./(1.-z) ) * m_Jmax;
}

double q_qg_IF::RejectionWeight(const double z,const double y,
				const double eta, const double scale) {
  return ( (1.-z)/(1.-z+y) - (1.-z*z)/2. ) * J(y,eta,scale)/m_Jmax;
}

double q_qg_IF::Z() {
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran.Get());
}

double q_qg_IF::J(const double z,const double eta,const double scale) { 
  
  p_pdf->Calculate(eta/z,scale);
  double fresh = p_pdf->GetXPDF(m_flavs[0]);
  
  p_pdf->Calculate(eta,scale);
  double old = p_pdf->GetXPDF(m_flavs[0]);

  return 1./z * fresh/old;
}

//###########################################################################
// q->gq
//###########################################################################
q_gq_IF::q_gq_IF(const ATOOLS::Flavour & flav,PDF::PDF_Base * pdf) : 
  m_Jmax(0.5), p_pdf(pdf)
{
  m_type     = cstp::IF;
  m_flavs[0] = flav;
  m_flavs[1] = Flavour(kf::gluon);
  m_flavs[2] = flav;
}

double q_gq_IF::operator() (const double z,const double y,
			    const double eta, const double scale) {
  return s_CF* ( 2./z - 2. +z ) * J(z,eta,scale);
}

double q_gq_IF::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  //if (scale<m_flavs[0].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1./xbj; 
  m_lastint = 2.*s_CF*log(zmax/zmin) * m_Jmax;
  return m_lastint;
}

double q_gq_IF::Overestimated(const double z,const double y) {
  return s_CF* ( 2./z ) * m_Jmax;
}

double q_gq_IF::RejectionWeight(const double z,const double y,
				const double eta, const double scale) {
  return ( 1. - z + z*z/2. ) * J(z,eta,scale)/m_Jmax;
}

double q_gq_IF::Z() {
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double q_gq_IF::J(const double z,const double eta,const double scale) { 
  
  p_pdf->Calculate(eta/z,scale);
  double fresh = p_pdf->GetXPDF(m_flavs[0]);
  
  p_pdf->Calculate(eta,scale);
  double old = p_pdf->GetXPDF(m_flavs[0]);

  return 1./z * fresh/old;
}

 
//###########################################################################
// g->gg
//###########################################################################
g_gg_IF::g_gg_IF(PDF::PDF_Base * pdf) : m_Jmax(0.5), p_pdf(pdf)
{
  m_type     = cstp::IF;
  m_flavs[0] = Flavour(kf::gluon);
  m_flavs[1] = Flavour(kf::gluon);
  m_flavs[2] = Flavour(kf::gluon);
}

double g_gg_IF::operator() (const double z,const double y,
			    const double eta, const double scale) {
  return s_CA * ( 1./(1.-z+y) + 1./z -2. +z*(1.-z) ) * J(z,eta,scale);
}

double g_gg_IF::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1./xbj; 
  m_lastint = s_CA * log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax))) * m_Jmax;
  return m_lastint;
}

double g_gg_IF::Overestimated(const double z,const double y) {
  return s_CA * ( 1./(1.-z) + 1./z ) * m_Jmax;
}

double g_gg_IF::RejectionWeight(const double z,const double y,
				const double eta, const double scale) {
  return z*(1.-z) * ( 1./(1.-z+z*y) + 1./z -2. + z*(1.-z) ) * J(z,eta,scale)/m_Jmax;
}

double g_gg_IF::Z() {
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double g_gg_IF::J(const double z,const double eta,const double scale) { 
  
  p_pdf->Calculate(eta/z,scale);
  double fresh = p_pdf->GetXPDF(m_flavs[0]);
  
  p_pdf->Calculate(eta,scale);
  double old = p_pdf->GetXPDF(m_flavs[0]);

  return 1./z * fresh/old;
}


//###########################################################################
// g->qq
//###########################################################################
g_qq_IF::g_qq_IF(const ATOOLS::Flavour & flav,PDF::PDF_Base * pdf) : 
  m_Jmax(0.5), p_pdf(pdf) 
{
  m_type     = cstp::IF;
  m_flavs[0] = Flavour(kf::gluon);
  m_flavs[1] = flav;
  m_flavs[2] = flav.Bar();
}

double g_qq_IF::operator() (const double z,const double y,
			    const double eta, const double scale) {
  return s_TR * (1.-2.*z*(1.-z)) * J(z,eta,scale);
}

double g_qq_IF::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = m_zmin = m_zmax = 0.; 
  if (scale<m_flavs[1].PSMass()) return m_lastint;
  m_zmin = zmin; m_zmax = zmax;
  m_Jmax = 1./xbj; 
  m_lastint = s_TR*(m_zmax-m_zmin) * m_Jmax;                                             
  return m_lastint;
}

double g_qq_IF::Overestimated(const double z,const double y) {
  return s_TR * m_Jmax;
}

double g_qq_IF::RejectionWeight(const double z,const double y,
				const double eta, const double scale) {
  return (1.-2.*z*(1.-z)) * J(z,eta,scale)/m_Jmax;
}

double g_qq_IF::Z() {
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double g_qq_IF::J(const double z,const double eta,const double scale) { 
  
  p_pdf->Calculate(eta/z,scale);
  double fresh = p_pdf->GetXPDF(m_flavs[0]);
  
  p_pdf->Calculate(eta,scale);
  double old = p_pdf->GetXPDF(m_flavs[0]);

  return 1./z * fresh/old;
}
