#include "SHRiMPS/Eikonals/Eikonal_Weights.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Histogram.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;


Eikonal_Weights::Eikonal_Weights() :
  m_lambda(MBpars.GetEikonalParameters().lambda),
  m_Delta(MBpars.GetEikonalParameters().Delta),
  m_originalY(MBpars.GetEikonalParameters().originalY),
  m_Ymax(MBpars.GetEikonalParameters().Ymax),
  m_bmax(MBpars.GetEikonalParameters().bmax),
  m_density(Rapidity_Density(m_Delta,m_lambda))
{ }

void Eikonal_Weights::SetEikonal(Omega_ik * eikonal) {
  m_density.SetEikonal(eikonal);
}

void Eikonal_Weights::
SetImpactParameters(const double & b1, const double & b2) {
  m_b1 = b1; m_b2 = b2;
  m_density.SetImpactParameters(m_b1,m_b2);
}

void Eikonal_Weights::
AddRapidities(Ladder * ladder,const double & ymin,const double & ymax) {
  size_t ngluons(m_density.NGluons(ymin,ymax));
  for (size_t i=0;i<ngluons;i++)
    ladder->AddRapidity(m_density.SelectRapidity(ymin,ymax));
}

colour_type::code Eikonal_Weights::
PropColour(const double & y1,const double & y2) {
  double singletwt(m_density.SingletWeight(y1,y2));
  double octetwt(m_density.OctetWeight(y1,y2));
  double tot(singletwt+octetwt);
  if (ran->Get()*tot>octetwt) return colour_type::singlet;
  return colour_type::octet;
}

double Eikonal_Weights::
WeightSingletOverOctet(const double & y1,const double & y2) {
  double singletwt(m_density.SingletWeight(y1,y2));
  double octetwt(m_density.OctetWeight(y1,y2));
  return singletwt/octetwt;
}   


void Eikonal_Weights::Test(const std::string & dirname) {
  double ymax(MBpars.GetEikonalParameters().Ymax);
  m_density.SetLambdaForTest(0.25);
  Histogram histo_ngluon(0,0.,20.,20);
  Histogram histo_ygluon(0,-10.,10.,100);
  long int ng(0), steps(100000);
  for (int i=0;i<steps;i++) {
    ng += m_density.NGluons(-ymax,ymax);
    histo_ngluon.Insert(double(ng)+0.5);
    histo_ygluon.Insert(m_density.SelectRapidity(-ymax,ymax)+0.1);
  }
  msg_Info()<<METHOD<<" yields <ng> = "<<(double(ng)/double(steps))
	    <<" vs. "<<m_density.MeanNGluons()<<".\n";
  string name_ngluon(string("NGluons_In_Ladder.dat"));
  histo_ngluon.Finalize();
  histo_ngluon.Output(dirname+"/"+name_ngluon);
  string name_ygluon(string("YGluons_In_Ladder.dat"));
  histo_ygluon.Finalize();
  histo_ygluon.Output(dirname+"/"+name_ygluon);
  Histogram histo_ydef(0,-10.,10.,100);
  double y(-10.+0.1);
  while (y<10.) {
    histo_ydef.Insert(y,m_density(y));
    y += 0.1;
  }
  string name_ydef(string("YGluons_In_Ladder_Direct.dat"));
  histo_ydef.Finalize();
  histo_ydef.Output(dirname+"/"+name_ydef);
}


// double Eikonal_Weights::
// MaximalEmissionProbability(const double & b1,const double & b2) 
// {
//   return m_Delta;
// }

// double Eikonal_Weights::
// EmissionWeight(const double & b1,const double & b2,const double & y,
// 	       const double & sup) {
//   if (y<-m_originalY || y>m_originalY || b1>m_bmax || b2>m_bmax) return 0.;
//   if (dabs(y)>m_Ymax) return 1.;
//   double term1 = ATOOLS::Max(1.e-12,m_lambda/2.*sup*(*p_Omegaik)(b1,b2,y));
//   double term2 = ATOOLS::Max(1.e-12,m_lambda/2.*sup*(*p_Omegaki)(b1,b2,y));
//   double absorption(1.);
//   switch (m_absorp) {
//   case absorption::factorial:
//     if (!ATOOLS::IsZero(term1) && !ATOOLS::IsZero(term2)) 
//       absorption = (1.-exp(-term1))/term1 * (1.-exp(-term2))/term2;
//     break;
//   case absorption::exponential:
//   default:
//     absorption = exp(-(term1+term2));
//     break;
//   }
//   return absorption;
// }

// double Eikonal_Weights::
// SingletWeight(const double & b1,const double & b2,
// 	      const double & y1,const double & y2,
// 	      const double & sup,const int & nbeam) {
//   double term   = m_singletwt*DeltaOmega(b1,b2,y1,y2,sup,nbeam); 
//   double weight = sqr(1.-exp(-term/2.));
//   return weight;
// }

// double Eikonal_Weights::
// OctetWeight(const double & b1,const double & b2,
// 	    const double & y1,const double & y2,
// 	    const double & sup,const int & nbeam) {
//   double term   = DeltaOmega(b1,b2,y1,y2,sup,nbeam); 
//   double weight = 1.-exp(-term);
//   return weight;
// }

// double Eikonal_Weights::
// RescatterProbability(const double & b1,const double & b2,
// 		     const double & y1,const double & y2,
// 		     const double & sup,const int & nbeam) {
//   double term   = DeltaOmega(b1,b2,y1,y2,sup,nbeam); 
//   double weight = 1.-exp(-term);
//   return weight;
// }

// double Eikonal_Weights::
// DeltaOmega(const double & b1,const double & b2,
// 	   const double & y1,const double & y2,
// 	   const double & sup,const int & nbeam) {
//   if (b1<0. || b1>m_bmax || b2<0. || b2>m_bmax)     return 0.;
//   if (dabs(y1)>m_originalY || dabs(y2)>m_originalY) return 0.;
//   double meany((y1+y2)/2.), ommaj, ommin;
//   if (meany<0.) {
//     ommaj = (y1<y2)?(*p_Omegaik)(b1,b2,y2):(*p_Omegaik)(b1,b2,y1);
//     ommin = (y1<y2)?(*p_Omegaik)(b1,b2,y1):(*p_Omegaik)(b1,b2,y2);
//   }
//   else {
//     ommaj = (y1<y2)?(*p_Omegaki)(b1,b2,y1):(*p_Omegaki)(b1,b2,y2);
//     ommin = (y1<y2)?(*p_Omegaki)(b1,b2,y2):(*p_Omegaki)(b1,b2,y1);
//   }
//   return sup*pow(m_lambda,2-nbeam)*dabs(ommaj-ommin)/(ommin);
// }

// double Eikonal_Weights::
// EffectiveIntercept(double b1,double b2,const double & y)
// {
//   if (b1<0. || b1>m_bmax || b2<0. || b2>m_bmax || 
//       dabs(y)>m_originalY) {
//     return 0.;
//   }
//   double res(m_Delta * exp(-m_lambda *
// 			   ((*p_Omegaik)(b1,b2,y)+(*p_Omegaki)(b1,b2,y))/2.));
//   return res;;
// }

// double Eikonal_Weights::
// Sum(const double & b1,const double & b2,const double & y){
//   if (dabs(y)>m_originalY) return 0.;
//   if (dabs(y)>m_Ymax) return 1.;
//   double term1 = (*p_Omegaik)(b1,b2,y)/p_ff1->FourierTransform(b1);
//   double term2 = (*p_Omegaki)(b1,b2,y)/p_ff2->FourierTransform(b2);

//   return term1+term2;
// }
