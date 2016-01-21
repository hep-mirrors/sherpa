#include "SHRiMPS/Event_Generation/Parton_Luminosity.H"
#include "SHRiMPS/Event_Generation/Inelastic_Event_Generator.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Parton_Luminosity::Parton_Luminosity(Beam_Remnant_Handler * beams) : 
  m_adjust(true), 
  m_Emin(2.), m_Ecms(rpa->gen.Ecms()), 
  m_smin(4.*sqr(m_Emin)), m_smintest(m_smin), m_smax(sqr(m_Ecms)),
  m_shatsteps(1000),
  m_kernel(m_smin,m_smax)
{
  for (size_t beam=0;beam<2;beam++) p_pdf[beam] = beams->GetPDF(beam);
  m_kernel.SetISR(p_pdf);
}

Parton_Luminosity::~Parton_Luminosity() {}

void Parton_Luminosity::FillGrids(Inelastic_Event_Generator * generator) {
  std::list<Omega_ik *> * eikonals(MBpars.GetEikonals());  
  for (std::list<Omega_ik *>::iterator eikiter=eikonals->begin();
       eikiter!=eikonals->end();eikiter++) {
    m_kernel.SetExponent(1.+(*eikiter)->EffectiveIntercept());
    m_smin    = m_smintest;
    m_kernel.SetSmin(m_smin);
    double maxdl, totallumi = CalculateTotalXSec(maxdl);
    if (m_adjust) {
      double s0ratio = pow(totallumi*rpa->Picobarn()/
			   generator->XSec(*eikiter),
			   1./(1.+(*eikiter)->EffectiveIntercept()));
      m_smin *= s0ratio;
      m_kernel.SetSmin(m_smin);
      totallumi = CalculateTotalXSec(maxdl);
    }
    m_maxdls[(*eikiter)]    = maxdl;
    m_deltalogs[(*eikiter)] = m_deltalog;
    m_smins[(*eikiter)]     = m_smin;
  }
}

double Parton_Luminosity::
CalculateTotalXSec(double & maxdl) {
  Gauss_Integrator integrator((&m_kernel));
  double intl(0.);
  m_deltalog = log(m_smax/m_smin-1.)/double(m_shatsteps+1);
  maxdl = 0.;
  double shat(0.),shat1(0.),ymax(0.),dl(0.),dl1(0.);
  for (int i=0;i<=m_shatsteps;i++) {
    shat = m_smin*exp(m_deltalog*i);
    m_kernel.Reset();
    m_kernel.SetShat(shat);
    ymax = -log(shat/m_smax)/2.;
    dl   = integrator.Integrate(-ymax,ymax,0.01,1)/shat;
    if (maxdl<m_kernel.MaxDL()) maxdl = m_kernel.MaxDL();
    if (i>0) intl += (shat-shat1)/(2.*shat)*(dl+dl1)/2.;
    shat1 = shat;
    dl1   = dl;
  }
  return intl;
}

double Parton_Luminosity::Kernel::operator()(double y) {
  double pref(sqrt(m_shat/m_smax)),x0(pref*exp(-y)),x1(pref*exp(y));
  double pdfs(0.);
  if (x0>p_pdf[0]->XMin() && x1>p_pdf[1]->XMin() && 
      x0<p_pdf[0]->XMax() && x1<p_pdf[1]->XMax()) {
    pdfs = p_pdf[0]->AllPartons(x0,0.)*p_pdf[1]->AllPartons(x1,0.);
  }
  double val(pdfs*Weight(m_shat));
  if (val>m_maxdl) m_maxdl = val;
  return val;
}

void Parton_Luminosity::SetEikonal(Omega_ik * eikonal) {
  p_eikonal  = eikonal;
  m_maxdl    = m_maxdls[p_eikonal];
  m_deltalog = m_deltalogs[p_eikonal];
  m_smin     = m_smins[p_eikonal];
  m_kernel.SetSmin(m_smin);
  m_eta      = p_eikonal->EffectiveIntercept();
}
