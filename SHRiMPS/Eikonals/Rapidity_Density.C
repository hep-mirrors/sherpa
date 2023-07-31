#include "SHRiMPS/Eikonals/Rapidity_Density.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;


//////////////////////////////////////////////////////////
//
// operator:
//       d<N>/dy = Delta * exp[-lambda/2 * ( Omega_{ik}(b_1,b_2,y) + Omega_{ki}(b_1,b_2,y) ) ]
//
//  Integrate:
//       <N> = int dy d<N>/dy
//
//  NGluons:
//       N = Poissonian(<N>) 
//
//  SelectRapidities:
//       selected according to d<N>/dy
//
//  DeltaOmega: 
//       delta Omega(y_1,y_2) = Omega_{ik}(b_1,b_2,y_2) / Omega_{ik}(b_1,b_2,y_1) - 1
//             assuming that y_2 > y_1   and (y_1+y_2) / 2 < 0.
//             otherwise - i.e. for (y_1+y_2)/2 > 0 we use Omega_{ki}
//
//  SingletWeight:
//       W_1(b_1,b_2,y_1,y_2) = [1-exp(-delta Omega/2)]^2
//
//  OctetWeight:
//       W_8(b_1,b_2,y_1,y_2) = 1-exp(-delta Omega)
//
/////////////////////////////////////////////////////////

Rapidity_Density::Rapidity_Density(const double & Delta,const double & lambda,
				   const double & Ymax,
				   const absorption::code & absorp) :
  m_Delta(Delta), m_lambda(lambda), m_Ymax(Ymax), m_b1(0.), m_b2(0.),
  m_absorp(absorp)
{}

void Rapidity_Density::Test(Omega_ik * eikonal) {
  SetEikonal(eikonal);
}

void Rapidity_Density::SetEikonal(Omega_ik * eikonal) {
  p_omegaik = eikonal->GetSingleTerm(0);
  p_omegaki = eikonal->GetSingleTerm(1);
}

void Rapidity_Density::
SetImpactParameters(const double & b1, const double & b2) {
  m_b1 = b1; m_b2 = b2; m_max = 0.; m_mean = 0.;
}

double Rapidity_Density::AbsorptionWeight(double y) {
  double O_ik = m_lambda/2.*(*p_omegaik)(m_b1,m_b2,y);
  double O_ki = m_lambda/2.*(*p_omegaki)(m_b1,m_b2,y);
  switch (m_absorp) {
  case absorption::truncatedExp:
    return 1.-(O_ik+O_ki);
    break;
  case absorption::exponential:
    return exp(-(O_ik+O_ki));
    break;
  case absorption::factorial:
  default:
    return (1.-exp(-O_ik))/O_ik * (1.-exp(-O_ki))/O_ki;
    break;
  }
  return 0.;
}

double Rapidity_Density::operator()(double y) {
  //return m_Delta;
  double result = m_Delta * AbsorptionWeight(y);
  if (result>m_max) m_max=result;
  return result;
}

double Rapidity_Density::Integrate(const double & ymin,const double & ymax) {
  ATOOLS::Gauss_Integrator integrator(this);
  return integrator.Integrate(ymin,ymax,1.e-5,1);
}

size_t Rapidity_Density::NGluons(const double & ymin,const double & ymax,const bool & rescatter) {
  if (rescatter) return ran->Poissonian(m_Delta * dabs(ymax-ymin));
  m_mean = dabs(Integrate(ymin,ymax));
  return ran->Poissonian(m_mean);
}

double Rapidity_Density::SelectRapidity(const double & ymin,const double & ymax) {
  double y;
  do {
    y = ymin+ran->Get()*(ymax-ymin);
  } while ((*this)(y)<MaxWeight()*ran->Get());
  return y;
}

double Rapidity_Density::RescatterProbability(const double & y1,const double & y2) {
  return 1.-exp(-DeltaOmega(y1,y2));
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
  double expo = int(dabs(y1)<m_Ymax)+int(dabs(y2)<m_Ymax);
  return pow(m_lambda,expo) * dabs(ommaj-ommin)/(ommin);
}

double Rapidity_Density::SingletWeight(const double & y1,const double & y2) {
  return sqr(1.-exp(-DeltaOmega(y1,y2)/2.));
}

double Rapidity_Density::OctetWeight(const double & y1,const double & y2) {
  return 1.-exp(-DeltaOmega(y1,y2));
}

double Rapidity_Density::EffectiveIntercept(const double & b1, const double & b2,const double & y) {
  if (dabs(y)>m_Ymax) return 0.;
  // watch below: there may be a factor 1/2 missing in the exponential
  return m_Delta * exp(-m_lambda * ( (*p_omegaik)(b1,b2,y)+(*p_omegaki)(b1,b2,y)) );
}


