#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "SHRiMPS/Beam_Remnants/Remnant_Handler.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;
  
Sigma_Partonic::Sigma_Partonic(const xs_mode::code & mode) :
  p_alphaS(NULL), m_mode(mode), m_fixflavour(true), 
  m_Ymax(MBpars.GetEikonalParameters().originalY),
  m_S(sqr(rpa->gen.Ecms())),
  m_eta(0.08), m_smin(1.), m_tmin(1.),
  m_accu(0.01), m_sigma(0.), m_maxdsigma(0.),
  m_Nmaxtrials(100), m_kinX_fails(0)
{ }

Sigma_Partonic::~Sigma_Partonic() {}

void Sigma_Partonic::Initialise(Remnant_Handler * remnants) {
  for (size_t beam=0;beam<2;beam++) p_pdf[beam]=remnants->GetPDF(beam);
  m_smin = Max(m_smin, m_S*p_pdf[0]->XMin()*p_pdf[1]->XMin());
  if (!Calculate()) {
    msg_Error()<<METHOD<<" fails: integration did not converge.  "
	       <<"Will exit the run.\n"
	       <<"   Encountered "<<m_kinX_fails<<" fails in creating good kinematics.\n";
    exit(1);
  }
}

const double Sigma_Partonic::MakeEvent() {
  for (size_t i=0;i<m_Nmaxtrials;i++) {
    m_dsigma = MakePoint() * dSigma();
    if (m_dsigma>m_maxdsigma*ran->Get()) {
      m_phi = 2.*M_PI*ran->Get();
      SelectFlavours(m_fixflavour); 
      return m_shat;
    }
  }
  return -1.;
}

void Sigma_Partonic::SelectFlavours(const bool & fixflavour) {
  m_flavs[0] = m_flavs[1] = Flavour(kf_gluon);
  if (fixflavour) return;
  for (size_t i=0;i<2;i++) {
    double disc = 0.999999 * m_xpdf[i]*ran->Get();
    for (list<Flavour>::const_iterator flit=p_pdf[i]->GetFlavours().begin();
	 flit!=p_pdf[i]->GetFlavours().end();flit++) {
      if (p_pdf[i]->XPDF((*flit))<1.e-4) continue;
      disc -= p_pdf[i]->XPDF((*flit)) * ColourFactor((*flit));
      if (disc<=0.) {
	m_flavs[i] = (*flit);
	break;
      }
    }
  }
}

const double Sigma_Partonic::
PDF(const size_t beam,const double & x,const double & Q2,const Flavour & flav) {
  if (x<1.e-6) return 0.;
  p_pdf[beam]->Calculate(x,Q2);
  return p_pdf[beam]->XPDF(flav);  
}

const bool Sigma_Partonic::Calculate() {
  size_t iter = 0, Npoints = 10000;
  long int N  = 0;
  double dsigma, sigma, res = 0., res2 = 0., accu = 1.e99;
  do {
    for (size_t i=0;i<Npoints;i++) {
      N++;
      dsigma = MakePoint() * dSigma();
      res   += dsigma;
      res2  += sqr(dsigma);
      if (dsigma>m_maxdsigma) m_maxdsigma = dsigma;
    }
    sigma = res/double(N);
    accu  = sqrt((res2/double(N) - sqr(sigma))/double(N))/sigma;
    if (accu<m_accu) {
      m_Nmaxtrials = Max(1000,int(10.*m_maxdsigma/sigma));
      msg_Out()<<METHOD<<" succeeds after "<<N<<" points:\n"
	       <<"  sigma = "<<(sigma*rpa->Picobarn()*1.e-9)<<" mb "
	       <<"+/- "<<(100.*accu)<<" %, "
	       <<"max value = "<<m_maxdsigma<<";\n"
	       <<"  expected unweighting efficiency = "
	       <<1./double(m_Nmaxtrials)<<" "
	       <<"from "<<sigma<<" and "<<m_maxdsigma<<" ==> "
	       <<m_Nmaxtrials<<"\n";
      return true;
    }
    iter++;
  } while (iter<100 && accu>m_accu);
  msg_Out()<<METHOD<<" integration after "<<N<<" points dos not converge:\n"
	   <<"   sigma = "<<(sigma*rpa->Picobarn()*1.e-9)<<" mb "
	   <<"+/- "<<(100.*accu)<<" %, "
	   <<"max value = "<<m_maxdsigma<<".\n"
	   <<"   Encountered "<<(double(m_kinX_fails)/double(100*Npoints)*100)
	   <<"% fails in creating good kinematics.\n";
  return false;
}

const double Sigma_Partonic::MakePoint() {
  double cosh2, pt2min, pt2max, ymax, rand = ran->Get();
  switch (m_mode) {
  case xs_mode::perturbative:
    for (size_t i=0;i<2;i++) {
      m_y[i] = m_Ymax*(-1.+2.*ran->Get());
    }
    if (m_y[0]==m_y[1])     { m_y[0]+=1.e-6; m_y[1]-=1.e-6; }
    else if (m_y[0]<m_y[1]) { double ysave = m_y[1]; m_y[1] = m_y[0]; m_y[0] = ysave;}
    m_yhat   = (m_y[0]+m_y[1])/2.;
    m_dy     = (m_y[0]-m_y[1])/2.;
    m_coshdy = cosh(m_dy);
    cosh2    = sqr(m_coshdy);
    pt2min   = m_smin/(4.*cosh2);
    pt2max   = m_S/(4.*cosh2);
    //m_pt2    = 1./(rand/pt2max + (1.-rand)/pt2min);
    m_pt2    = (pt2min+m_tmin)*pow((pt2max+m_tmin)/(pt2min+m_tmin),rand)-pt2min;
    m_shat   = 4.*cosh2*m_pt2;
    m_that   = -m_pt2*(1.+exp(-m_dy));
    //return sqr(2.*m_Ymax) * (1./pt2min-1./pt2max) / (8.*M_PI*m_shat);
    return sqr(2.*m_Ymax) * log((pt2max+m_tmin)/(pt2min+m_tmin)) / (8.*M_PI*m_shat);
  case xs_mode::integrated:
  case xs_mode::Regge:
    m_shat   = m_smin * pow(m_S/m_smin,ran->Get());
    ymax     = 1./2.*log(m_S/m_shat);
    m_yhat   = ymax * (-1. + 2.*rand);
    return  log(m_S/m_smin) * (2.*ymax);
  default:
    break;
  }
  return 0.;
}

const double Sigma_Partonic::dSigma() {
  double flux = 1/(2.*m_shat), Enorm = sqrt(m_shat/m_S), scale = m_that;
  for (size_t i=0;i<2;i++) {
    m_x[i]    = Enorm * (i==0 ? exp(m_yhat) : exp(-m_yhat));
    if (m_x[i]<p_pdf[i]->XMin() || m_x[i]>0.9999999) {
      m_kinX_fails++;
      return 0.;
    }
    m_xpdf[i] = 0.;
    p_pdf[i]->Calculate(m_x[i],scale);
    for (list<Flavour>::const_iterator flit=p_pdf[i]->GetFlavours().begin();
	 flit!=p_pdf[i]->GetFlavours().end();flit++) {
      if (p_pdf[i]->XPDF((*flit))>1.e-4) 
	m_xpdf[i] += p_pdf[i]->XPDF((*flit)) * ColourFactor((*flit));
    }
  }
  return flux * m_xpdf[0]*m_xpdf[1] * ME2(m_shat,m_that,scale);
}

const double Sigma_Partonic::dSigma(const double & shat,const double & yhat) {
  double flux = 1/(2.*m_shat), Enorm = sqrt(m_shat/m_S), scale = 0.;
  for (size_t i=0;i<2;i++) {
    m_x[i]    = Enorm * (i==0 ? exp(m_yhat) : exp(-m_yhat));
    if (m_x[i]<1.e-6 || m_x[i]>0.9999999) return 0.;
    m_xpdf[i] = 0.;
    p_pdf[i]->Calculate(m_x[i],scale);
    for (list<Flavour>::const_iterator flit=p_pdf[i]->GetFlavours().begin();
	 flit!=p_pdf[i]->GetFlavours().end();flit++) {
      if (p_pdf[i]->XPDF((*flit))<1.e-4) continue;
      m_xpdf[i] += p_pdf[i]->XPDF((*flit)) * ColourFactor((*flit));
    }
  }
  return flux * m_xpdf[0]*m_xpdf[1] * ME2(m_shat,m_that,scale);
}

const double Sigma_Partonic::ME2(const double & shat,const double & that,const double & scale) {
  switch (m_mode) {
  case xs_mode::perturbative:
    return sqr(4.*M_PI*(*p_alphaS)(scale))*(sqr(shat)+sqr(shat+that))/sqr(-that+m_tmin);
  case xs_mode::integrated:
    return 4.*M_PI*sqr((*p_alphaS)(scale))*sqr(shat)/(m_tmin*(shat+m_tmin));
  case xs_mode::Regge:
    return pow( (shat/m_smin), 1.+m_eta);
  default:
    break;
  }
  return 0.;
}
