#include "SHRiMPS/Eikonals/Rapidity_Density.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;


Rapidity_Density::Rapidity_Density(const double & Delta,const double & lambda) :
  m_Delta(Delta), m_lambda(lambda)
{}

void Rapidity_Density::SetEikonal(Omega_ik * eikonal) {
  p_omegaik = eikonal->GetSingleTerm(0);
  p_omegaki = eikonal->GetSingleTerm(1);
}

void Rapidity_Density::
SetImpactParameters(const double & b1, const double & b2) {
  m_b1 = b1; m_b2 = b2; m_max = 0.; m_mean = 0.;
}

double Rapidity_Density::operator()(double y) {
  double result = m_Delta;
  //double arg    = (*p_omegaik)(m_b1,m_b2,y)+(*p_omegaki)(m_b1,m_b2,y);
  //double result = m_Delta * exp(-m_lambda/2.*arg);
  if (result>m_max) m_max=result;
  return result;
}

double Rapidity_Density::Integrate(const double & ymin,const double & ymax) {
  ATOOLS::Gauss_Integrator integrator(this);
  return integrator.Integrate(ymin,ymax,1.e-5,1);
}

size_t Rapidity_Density::NGluons(const double & ymin,const double & ymax) {
  m_mean = Integrate(ymin,ymax);
  return ran->Poissonian(m_mean);
}

double Rapidity_Density::
SelectRapidity(const double & ymin,const double & ymax) {
  double y;
  do {
    y = ymin+ran->Get()*(ymax-ymin);
  } while ((*this)(y)<MaxWeight()*ran->Get());
  return y;
}

double Rapidity_Density::MaxWeight() {
  return m_max;
}

double Rapidity_Density::DeltaOmega(const double & y1,const double & y2) {
  double meany((y1+y2)/2.), ommaj, ommin;
  if (meany<0.) {
    ommaj = (y1<y2)?(*p_omegaik)(m_b1,m_b2,y2):(*p_omegaik)(m_b1,m_b2,y1);
    ommin = (y1<y2)?(*p_omegaik)(m_b1,m_b2,y1):(*p_omegaik)(m_b1,m_b2,y2);
  }
  else {
    ommaj = (y1<y2)?(*p_omegaki)(m_b1,m_b2,y1):(*p_omegaki)(m_b1,m_b2,y2);
    ommin = (y1<y2)?(*p_omegaki)(m_b1,m_b2,y2):(*p_omegaki)(m_b1,m_b2,y1);
  }
  return /*sqr(m_lambda)*/dabs(ommaj-ommin)/(ommin);
}

double Rapidity_Density::SingletWeight(const double & y1,const double & y2) {
  return sqr(1.-exp(-DeltaOmega(y1,y2)/2.));
}

double Rapidity_Density::OctetWeight(const double & y1,const double & y2) {
  return 1.-exp(-DeltaOmega(y1,y2));
}



