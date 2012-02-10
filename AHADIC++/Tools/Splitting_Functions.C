#include "AHADIC++/Tools/Splitting_Functions.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;

ZForm::code AHADIC::DefineZForm(const int & zform) {
  switch (zform) {
  case 2:
    return ZForm::splitkernel_flat;
  case 1:
    return ZForm::splitkernel;
  case 0:
  default:
    break;
  }
  return ZForm::flat;
}


Splitting_Functions::
Splitting_Functions(const ZForm::code & zform,const int & masstreatment) :
  m_zform(zform), m_masstreatment(masstreatment),
  m_alpha(2. /*hadpars->Get(std::string(""))*/), 
  m_sigma(hadpars->Get(std::string("pt02"))), 
  m_kappa(hadpars->Get(std::string("P_qg_Exponent")))
{
}

double Splitting_Functions::
SelectZ(const double & zmin,const double & zmax,const double & expo,
	const bool & glusplit,const bool & leading) {
  switch (m_zform) {
  case ZForm::splitkernel_flat:
    if (!leading) break;
  case ZForm::splitkernel:
    return (glusplit?
	    SelectZFromG2QQSplittingFunction(zmin,zmax):
	    SelectZFromQ2QGSplittingFunction(zmin,zmax,expo));
  case ZForm::flat:
  default:
    break;
  }
  return zmin + ran->Get()*(zmax-zmin);
}

double Splitting_Functions::
SelectZFromG2QQSplittingFunction(const double & zmin,const double & zmax) {  
  double z;
  do {
    z = zmin+(zmax-zmin)*ran->Get();
  } while (4.*z*(1.-z)<ran->Get());
  return z;
}

double Splitting_Functions::
SelectZFromQ2QGSplittingFunction(const double & zmin,const double & zmax,
				 const double & expo) {
  double z,rn;
  do {
    rn = ran->Get();
    if (expo==1.) 
      z = (1.-(1.-zmin)*pow((1.-zmax)/(1.-zmin),rn));
    else
      z = 1.-pow(rn*pow(1.-zmax,1.-expo)+
		 (1.-rn)*pow(1.-zmin,1.-expo),1./(1.-expo));
  } while (false || 1.+z*z<2.*ran->Get());
  return z;
}

double Splitting_Functions::
Weight(const double & y,const double & z,const double & Q2,
       const double & m12, const double & m22, const double & m32,
       const bool & glusplit,const bool & leading) {
  double weight = J(y,Q2,m12,m22,m32,glusplit);
  switch (m_zform) {
  case ZForm::splitkernel_flat:
    if (!leading) break; 
  case ZForm::splitkernel: 
    {
      double Qt2(Q2-m12-m22-m32),pipj(y*Qt2/2.);
      double vijk(sqrt(sqr(2.*m12+Qt2*(1.-y))-4.*m12)/(Qt2*(1.-y)));
      if (glusplit) {
	double viji(sqrt(sqr(Qt2*y)-4.*m22*m32)/(Qt2*y+2.*m32));
	double frac((2.*m32+Qt2*y)/(2.*(m32+m22+Qt2*m_y)));
	double zm(frac*(1.-viji*vijk)),zp(frac*(1.+viji*vijk));
	weight *= 
	  (1.-2.*(z*(1.-z)-zp*zm))/vijk * 
	  sqr(Qt2)/(Qt2+(m22+m32)/y)/sqrt(lambda(Q2,0,m12));
      }
      else {
	double vtijk(sqrt(lambda(Q2,m32,m12))/(Q2-m32-m12));
	weight *= 
	  pow((1.-z)*Max(0.,1./(1.-z+z*y) - 
			 vtijk/vijk *(1.+z+m32/pipj)/2.),m_kappa);
      }
    }
    break;
  case ZForm::flat:
  default:
    break;
  }


  return weight;
}

double Splitting_Functions::
J(const double & y,const double & Q2,
  const double & m12, const double & m22, const double & m32,
  const bool & glusplit) {
  double term(sqr(Q2-m12-m32)), mij2(m32);
  if (glusplit) mij2 = 0.;
  return (1.-y)*term/((Q2-m12-m32+(m32-mij2)/y)*sqrt(lambda(Q2,mij2,m12)));
}

double Splitting_Functions::
lambda(const double & x,const double & y,const double & z) {
  return x*x+y*y+z*z-2.*(x*y+y*z+z*z);
}

double Splitting_Functions::Integrated(const double & zmin,const double & zmax,
				       const bool & glusplit) {
  if (!glusplit) {
    switch (m_zform) {
    case ZForm::splitkernel_flat:
    case ZForm::splitkernel:
      return IntegratedFromQ2QGSplittingFunction(zmin,zmax);
    case ZForm::flat:
    default:
      break;
    }
  }
  return zmax-zmin;
}

double Splitting_Functions::
IntegratedFromQ2QGSplittingFunction(const double & zmin,const double & zmax) {
  double integrated(zmax-zmin);
  if (m_kappa==1.) 
    integrated = log((1.-zmin)/(1.-zmax));
  else 
    integrated = (pow(1.-zmin,1.-m_kappa)-pow(1.-zmax,1.-m_kappa))/(1.-m_kappa);
  return integrated;
}















