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
  
Sigma_Partonic::Sigma_Partonic(MODEL::Strong_Coupling * alphaS,
			       const xs_mode::code & mode) :
  p_alphaS(alphaS), m_mode(mode), m_fixflavour(true), 
  m_Ymax(MBpars.GetEikonalParameters().originalY),
  m_S(sqr(rpa->gen.Ecms())),
  m_eta(0.08), m_smin(1.), m_tmin(0.),
  m_accu(0.005), m_sigma(0.), m_maxdsigma(0.),
  m_Nmaxtrials(100)
{}

Sigma_Partonic::~Sigma_Partonic() {}

void Sigma_Partonic::Initialise(Remnant_Handler * remnants) {
  for (size_t beam=0;beam<2;beam++) p_pdf[beam]=remnants->GetPDF(beam);
  if (!Calculate()) {
    msg_Error()<<METHOD<<" fails: integration did not converge.  "
	       <<"Will exit the run.\n";
    exit(1);
  }
}

const double Sigma_Partonic::MakeEvent() {
  double Jac, dsigma;
  for (size_t i=0;i<m_Nmaxtrials;i++) {
    Jac      = MakePoint(m_shat,m_yhat); 
    m_dsigma = dSigma(m_shat,m_yhat);
    if (Jac * m_dsigma>m_maxdsigma*ran->Get()) {
      SelectFlavours(m_fixflavour); 
      return m_shat;
    }
  }
  return -1.;
}

void Sigma_Partonic::SelectFlavours(const bool & fixflavour) {
  //msg_Out()<<"         "<<METHOD<<"("<<fixflavour<<")\n";
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
  double s, y, Jac;
  do {
    for (size_t i=0;i<Npoints;i++) {
      N++;
      Jac    = MakePoint(s,y); 
      dsigma = Jac * dSigma(s,y);
      res   += dsigma;
      res2  += sqr(dsigma);
      if (dsigma>m_maxdsigma) m_maxdsigma = dsigma;
    }
    sigma = res/double(N);
    accu  = sqrt((res2/double(N) - sqr(sigma))/double(N))/sigma;
    if (accu<m_accu) {
      m_Nmaxtrials = int(1000.*sigma/m_maxdsigma);
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
	   <<"max value = "<<m_maxdsigma<<".\n";
  return false;
}

const double Sigma_Partonic::MakePoint(double & s,double & y) {
  s           = m_smin * pow(m_S/m_smin,ran->Get());
  double ymax = 1./2.*log(m_S/s);
  y           = ymax * (-1. + 2.*ran->Get());
  double Jac  = log(m_S/m_smin) * (2.*ymax);
  return Jac;
}

const double Sigma_Partonic::dSigma(const double & s,const double & y) {
  double flux = 1/(2.*s), Enorm = sqrt(s/m_S)/2., scale = 0.;
  for (size_t i=0;i<2;i++) {
    m_x[i]    = Enorm * (i==0 ? exp(y) : exp(-y));
    if (m_x[i]<1.e-6 || m_x[i]>0.9999999) return 0.;
    m_xpdf[i] = 0.;
    p_pdf[i]->Calculate(m_x[i],scale);
    for (list<Flavour>::const_iterator flit=p_pdf[i]->GetFlavours().begin();
	 flit!=p_pdf[i]->GetFlavours().end();flit++) {
      if (p_pdf[i]->XPDF((*flit))<1.e-4) continue;
      m_xpdf[i] += p_pdf[i]->XPDF((*flit)) * ColourFactor((*flit));
    }
  }
  return flux * m_xpdf[0]*m_xpdf[1] * ME2(s);
}

const double Sigma_Partonic::ME2(const double & s) {
  double me2 = 0.;
  switch (m_mode) {
  case xs_mode::perturbative:
  case xs_mode::Regge:
    me2 = pow(s/m_smin,1.+m_eta);
    break;
  default:
    break;
  }
  return me2;
}
