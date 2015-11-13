#include "SHRiMPS/Eikonals/Eikonal_Weights.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;


Eikonal_Weights::Eikonal_Weights(const Eikonal_Parameters & params) 
{
  m_lambda    = params.lambda;
  m_Delta     = params.Delta;
  m_originalY = params.originalY;
  m_Ymax      = params.Ymax;
  m_bmax      = params.bmax;
}

void Eikonal_Weights::SetSingleOmegaTerms(Eikonal_Contributor * omegaik,
					  Eikonal_Contributor * omegaki) {
  p_Omegaik = omegaik;
  p_Omegaki = omegaki;
}


double Eikonal_Weights::
MaximalEmissionProbability(const double & b1,const double & b2) 
{
  return m_Delta;
}

double Eikonal_Weights::
EmissionWeight(const double & b1,const double & b2,const double & y,
	       const double & sup) {
  if (y<-m_originalY || y>m_originalY || b1>m_bmax || b2>m_bmax) return 0.;
  if (dabs(y)>m_Ymax) return 1.;
  double term1 = ATOOLS::Max(1.e-12,m_lambda/2.*sup*(*p_Omegaik)(b1,b2,y));
  double term2 = ATOOLS::Max(1.e-12,m_lambda/2.*sup*(*p_Omegaki)(b1,b2,y));
  double absorption(1.);
  switch (m_absorp) {
  case absorption::factorial:
    if (!ATOOLS::IsZero(term1) && !ATOOLS::IsZero(term2)) 
      absorption = (1.-exp(-term1))/term1 * (1.-exp(-term2))/term2;
    break;
  case absorption::exponential:
  default:
    absorption = exp(-(term1+term2));
    break;
  }
  return absorption;
}

double Eikonal_Weights::
SingletWeight(const double & b1,const double & b2,
	      const double & y1,const double & y2,
	      const double & sup,const int & nbeam) {
  double term   = m_singletwt*DeltaOmega(b1,b2,y1,y2,sup,nbeam); 
  double weight = sqr(1.-exp(-term/2.));
  return weight;
}

double Eikonal_Weights::
OctetWeight(const double & b1,const double & b2,
	    const double & y1,const double & y2,
	    const double & sup,const int & nbeam) {
  double term   = DeltaOmega(b1,b2,y1,y2,sup,nbeam); 
  double weight = 1.-exp(-term);
  return weight;
}

double Eikonal_Weights::
RescatterProbability(const double & b1,const double & b2,
		     const double & y1,const double & y2,
		     const double & sup,const int & nbeam) {
  double term   = DeltaOmega(b1,b2,y1,y2,sup,nbeam); 
  double weight = 1.-exp(-term);
  return weight;
}

double Eikonal_Weights::
DeltaOmega(const double & b1,const double & b2,
	   const double & y1,const double & y2,
	   const double & sup,const int & nbeam) {
  if (b1<0. || b1>m_bmax || b2<0. || b2>m_bmax)     return 0.;
  if (dabs(y1)>m_originalY || dabs(y2)>m_originalY) return 0.;
  double meany((y1+y2)/2.), ommaj, ommin;
  if (meany<0.) {
    ommaj = (y1<y2)?(*p_Omegaik)(b1,b2,y2):(*p_Omegaik)(b1,b2,y1);
    ommin = (y1<y2)?(*p_Omegaik)(b1,b2,y1):(*p_Omegaik)(b1,b2,y2);
  }
  else {
    ommaj = (y1<y2)?(*p_Omegaki)(b1,b2,y1):(*p_Omegaki)(b1,b2,y2);
    ommin = (y1<y2)?(*p_Omegaki)(b1,b2,y2):(*p_Omegaki)(b1,b2,y1);
  }
  return sup*pow(m_lambda,2-nbeam)*dabs(ommaj-ommin)/(ommin);
}

double Eikonal_Weights::
EffectiveIntercept(double b1,double b2,const double & y)
{
  if (b1<0. || b1>m_bmax || b2<0. || b2>m_bmax || 
      dabs(y)>m_originalY) return 0.;
  return m_Delta * 
    exp(-m_lambda * ((*p_Omegaik)(b1,b2,y)+(*p_Omegaki)(b1,b2,y))/2.);
}

// double Eikonal_Weights::
// Sum(const double & b1,const double & b2,const double & y){
//   if (dabs(y)>m_originalY) return 0.;
//   if (dabs(y)>m_Ymax) return 1.;
//   double term1 = (*p_Omegaik)(b1,b2,y)/p_ff1->FourierTransform(b1);
//   double term2 = (*p_Omegaki)(b1,b2,y)/p_ff2->FourierTransform(b2);

//   return term1+term2;
// }
