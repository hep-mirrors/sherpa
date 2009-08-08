#include "AHADIC++/Tools/Splitting_Functions.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;

ZForm::code AHADIC::DefineZForm(const int & zform) {
  switch (zform) {
  case 30:
    return ZForm::fragmentation;
  case 3:
    return ZForm::splitkernel_frag;
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


Splitting_Functions::Splitting_Functions(const ZForm::code & zform,const int & masstreatment) :
  m_zform(zform), m_masstreatment(masstreatment),
  m_alpha(2. /*hadpars.Get(std::string(""))*/), 
  m_sigma(hadpars.Get(std::string("pt02"))), 
  m_kappa(hadpars.Get(std::string("P_qg_Exponent")))
{
}

double Splitting_Functions::SelectZ(const double & zmin,const double & zmax,
				    const bool & glusplit,const bool & leading) {
  if (glusplit) {
    switch (m_zform) {
    case ZForm::fragmentation:
    case ZForm::splitkernel_frag:
    case ZForm::splitkernel_flat:
    case ZForm::splitkernel:
      if (leading) return SelectZFromG2QQSplittingFunction(zmin,zmax);
      return zmin + ran.Get()*(zmax-zmin);
    case ZForm::flat:
    default:
      break;
    }
  }
  else {
    switch (m_zform) {
    case ZForm::fragmentation:
      return SelectZFromQ2QGFragmentationFunction(zmin,zmax);
    case ZForm::splitkernel_frag:
      if (leading) return SelectZFromQ2QGSplittingFunction(zmin,zmax);
      return SelectZFromQ2QGFragmentationFunction(zmin,zmax);
    case ZForm::splitkernel_flat:
      if (leading) return SelectZFromQ2QGSplittingFunction(zmin,zmax);
      return zmin + ran.Get()*(zmax-zmin);
    case ZForm::splitkernel:
      return SelectZFromQ2QGSplittingFunction(zmin,zmax);
    case ZForm::flat:
    default:
      break;
    }
  }
  return zmin + ran.Get()*(zmax-zmin);
}

double Splitting_Functions::
SelectZFromG2QQSplittingFunction(const double & zmin,const double & zmax) {  
  return zmin+(zmax-zmin)*ran.Get();
}

double Splitting_Functions::
SelectZFromQ2QGSplittingFunction(const double & zmin,const double & zmax) {
  if (m_kappa==1.) {
    double z(1.-(1.-zmin)*pow((1.-zmax)/(1.-zmin),ran.Get()));
    //std::cout<<METHOD<<" yields z = "<<z<<" in ["<<zmin<<", "<<zmax<<"]."<<std::endl;
    return z;
  }
  double rn = ran.Get();
  return 1.-pow(rn*pow(1.-zmax,1.-m_kappa)+(1.-rn)*pow(1.-zmin,1.-m_kappa),1./(1.-m_kappa));
}

double Splitting_Functions::
SelectZFromQ2QGFragmentationFunction(const double & zmin,const double & zmax) {
  return 1.-(1.-zmin)*pow((1.-zmax)/(1.-zmin),ran.Get());
}


double Splitting_Functions::Weight(const double & scale2,const double & z,
				   const bool & glusplit,const bool & leading) {
  if (glusplit) {
    switch (m_zform) {
    case ZForm::fragmentation:
    case ZForm::splitkernel_frag:
    case ZForm::splitkernel_flat:
      if (leading) return J(z,glusplit)*WeightForG2QQSplittingFunction(scale2,z);
      break;
    case ZForm::splitkernel:
      return J(z,glusplit)*WeightForG2QQSplittingFunction(scale2,z);
    case ZForm::flat:
    default:
      break;
    }
  }
  else {
    switch (m_zform) {
    case ZForm::fragmentation:
      return J(z,glusplit)*WeightForQ2QGFragmentationFunction(scale2,z);
    case ZForm::splitkernel_frag:
      if (leading) return J(z,glusplit)*WeightForQ2QGSplittingFunction(scale2,z);
      else return J(z,glusplit)*WeightForQ2QGFragmentationFunction(scale2,z);
    case ZForm::splitkernel_flat:
      if (leading)  return J(z,glusplit)*WeightForQ2QGSplittingFunction(scale2,z);
      break;
    case ZForm::splitkernel:
      return J(z,glusplit)*WeightForQ2QGSplittingFunction(scale2,z);
    case ZForm::flat:
    default:
      break;
    }
  }
  return J(z,glusplit);
}

double Splitting_Functions::
WeightForG2QQSplittingFunction(const double & scale2,const double & z) {  
  double Qt2 = m_Q2-m_m1_2-m_m2_2-m_m3_2;
  double wt  = (1.-2.*(z*(1.-z)-m_zp*m_zm))/m_vijk;
  wt *= sqr(Qt2)/(Qt2+(m_m2_2+m_m3_2)/m_y)/sqrt(lambda(m_Q2,0,m_m1_2));
  return wt;
}


double Splitting_Functions::
WeightForQ2QGSplittingFunction(const double & scale2,const double & z) {
  double vtijk = sqrt(lambda(m_Q2,m_m3_2,m_m1_2))/(m_Q2-m_m3_2-m_m1_2);
  double wt    = Max(0.,1./(1.-z+z*m_y) - vtijk/m_vijk * (1.+z+m_m3_2/m_pipj)/2.);
  return pow((1.-z)*wt,m_kappa);
}

double Splitting_Functions::
WeightForQ2QGFragmentationFunction(const double & scale2,const double & z) {
  return pow(1.-z,m_kappa)*pow(z,m_alpha) * exp(-scale2/(m_sigma*(1.-z)));
}

double Splitting_Functions::J(const double & z,const bool & glusplit) {
  double term(sqr(m_Q2-m_m1_2-m_m3_2)), mij2(m_m3_2);
  if (glusplit) mij2 = 0.;
  return (1.-m_y)*term/((m_Q2-m_m1_2-m_m3_2+(m_m3_2-mij2)/m_y)*sqrt(lambda(m_Q2,mij2,m_m1_2)));				    
}

double Splitting_Functions::lambda(const double & x,const double & y,const double & z) {
  return x*x+y*y+z*z-2.*(x*y+y*z+z*z);
}

double Splitting_Functions::Integrated(const double & zmin,const double & zmax,
				       const bool & glusplit) {
  if (!glusplit) {
    switch (m_zform) {
    case ZForm::fragmentation:
    case ZForm::splitkernel_frag:
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















