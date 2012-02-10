#include "SINIC++/Cross_Sections/Parton_Luminosity.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SINIC;
using namespace ATOOLS;

Parton_Luminosity::Parton_Luminosity(const double & Emin,const double & Ecms,
				     const int & shatsteps,const int & test) : 
  m_adjust(true), 
  m_Emin(Emin), m_Ecms(Ecms), 
  m_smin(4.*sqr(m_Emin)), m_smintest(m_smin), 
  m_smax(sqr(m_Ecms)), m_shatsteps(shatsteps), 
  m_test(test)
{ }

Parton_Luminosity::~Parton_Luminosity() {}

void Parton_Luminosity::FillGrids(std::list<Omega_ik *> * eikonals) {
  msg_Tracking()<<METHOD<<" : fill luminosity grids."<<std::endl;
  m_kernel = Kernel(m_smin,m_smax);
  m_kernel.SetISR(p_pdf);

  double totallumi, maxdl, s0ratio;
  bool killedpdfs(false);

  for (std::list<Omega_ik *>::iterator eikiter=eikonals->begin();
       eikiter!=eikonals->end();eikiter++) {
    m_kernel.SetExponent(1.+(*eikiter)->EffectiveIntercept());
    m_smin    = m_smintest;
    m_kernel.SetSmin(m_smin);
    totallumi = CalculateTotalXSec(maxdl);
    if (m_adjust) {
      msg_Tracking()<<"   Total xsection for eikonal_["
		    <<(*eikiter)->FF1()->Number()<<", "
		    <<(*eikiter)->FF2()->Number()<<"] "
		    <<"with eta = "<<(*eikiter)->EffectiveIntercept()<<" : "
		    <<totallumi*rpa->Picobarn()/1.e9<<" mbarn ("
		    <<(*eikiter)->Sigma_Inelastic()/1.e9
		    <<" mbarn); "
		    <<"deltalog = "<<m_deltalog<<", smin = "
		    <<m_smin<<"."<<std::endl;
      s0ratio   = pow(totallumi*rpa->Picobarn()/
		      (*eikiter)->Sigma_Inelastic(),
		      1./(1.+(*eikiter)->EffectiveIntercept()));
      m_smin *= s0ratio;
      m_kernel.SetSmin(m_smin);
      msg_Tracking()<<"         --> multiply m_smin with "<<s0ratio<<" = "
		    <<"pow("<<(totallumi*rpa->Picobarn()/
			       (*eikiter)->Sigma_Inelastic())<<"/"
		    <<1./(1.+(*eikiter)->EffectiveIntercept())<<") = "<<m_smin
		    <<"."<<std::endl;
      totallumi = CalculateTotalXSec(maxdl);
      msg_Info()<<"   After adjustment, total inelastic cross section "
		<<"for Omega_{"<<(*eikiter)->FF1()->Number()
		<<(*eikiter)->FF2()->Number()<<"} yields "
		<<totallumi*rpa->Picobarn()/1.e9<<" mbarn ("
		<<(*eikiter)->Sigma_Inelastic()/1.e9
		<<" mbarn); smin = "<<m_smin<<"."<<std::endl;
    }
    else {
      msg_Info()<<"   Without adjustment, total inelastic cross section "
		<<"for Omega_{"<<(*eikiter)->FF1()->Number()
		<<(*eikiter)->FF1()->Number()<<"} yields "
		<<totallumi*rpa->Picobarn()/1.e9<<" mbarn ("
		<<(*eikiter)->Sigma_Inelastic()/1.e9
		<<" mbarn); smin = "<<m_smin<<"."<<std::endl;
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
