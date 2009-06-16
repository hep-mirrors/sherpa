#include "AMISIC++/Model/Profile_Function_Base.H"

#include "AMISIC++/Model/Profile_Function.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"

using namespace AMISIC;

std::ostream &AMISIC::operator<<(std::ostream &ostr,const pft::code code)
{
  switch (code) {
  case pft::none:            return ostr<<"None";
  case pft::flat:            return ostr<<"Flat";
  case pft::exponential:     return ostr<<"Exponential";
  case pft::gaussian:        return ostr<<"Gaussian";
  case pft::double_gaussian: return ostr<<"Double Gaussian";
  }
  return ostr;
}

Differential_Overlap::Differential_Overlap(Profile_Function_Base *const owner):
  p_owner(owner) {}

double Differential_Overlap::operator()(const double b)
{
  return 2.0*M_PI*b*p_owner->KFactor()*(*p_owner)(b);
}

Interaction_Probability::Interaction_Probability(Profile_Function_Base *const owner):
  p_owner(owner) {}

double Interaction_Probability::operator()(const double b)
{
  return 2.0*M_PI*b*(1.0-exp(-p_owner->KFactor()*(*p_owner)(b)));
}

Profile_Function_Base::Profile_Function_Base(const pft::code code,
					     const double bmin,const double bmax):
  p_overlap(new Differential_Overlap(this)),
  p_probability(new Interaction_Probability(this)),
  m_type(code), 
  m_bmin(bmin),
  m_bmax(bmax),
  m_omin(0.0),
  m_omax(0.0),
  m_kfactor(1.0),
  m_omean(1.0),
  m_norm(1.0) {}

Profile_Function_Base::~Profile_Function_Base()
{
  delete p_probability;
  delete p_overlap;
}

Profile_Function_Base *Profile_Function_Base::SelectProfile(const std::string &type,
							    const std::vector<double> &parameters)
{
  Profile_Function_Base *profile=NULL;
  if ((profile=CreateProfile<Double_Gaussian_Profile>(type,parameters))!=NULL);
  else if ((profile=CreateProfile<Gaussian_Profile>(type,parameters))!=NULL);
  else if ((profile=CreateProfile<Exponential_Profile>(type,parameters))!=NULL);
  else profile=CreateProfile<Flat_Profile>(type,parameters);
  return profile;
}

#ifdef DEBUG__Profile_Function_Base
class Test_Func_O: public ATOOLS::Function_Base {
private:
  Profile_Function_Base *p_owner;
public:
  Test_Func_O(Profile_Function_Base *profile):p_owner(profile) {}
  double operator()(const double b) 
  { 
    return p_owner->KFactor()*(*p_owner)(b)*(*p_owner->Probability())(b);
  }
};
class Test_Func_F: public ATOOLS::Function_Base {
private:
  Profile_Function_Base *p_owner;
public:
  Test_Func_F(Profile_Function_Base *profile):p_owner(profile) {}
  double operator()(const double b) 
  { 
    return (*p_owner)(b)/p_owner->OMean()*(*p_owner->Probability())(b);
  }
};
#endif

bool Profile_Function_Base::CalculateOMean(const double ratio)
{
  ATOOLS::Gauss_Integrator *gausso = new ATOOLS::Gauss_Integrator(p_overlap);
  ATOOLS::Gauss_Integrator *gaussp = new ATOOLS::Gauss_Integrator(p_probability);
  double k1=ratio, k2=ratio;
  m_kfactor=k2;
  double ratio2=gausso->Integrate(m_bmin,m_bmax,1.0e-5);
  ratio2/=gaussp->Integrate(m_bmin,m_bmax,1.0e-5);
  m_kfactor=2.0*k2;
  do {
    double ratio1=ratio2;
    ratio2=gausso->Integrate(m_bmin,m_bmax,1.0e-5);
    ratio2/=gaussp->Integrate(m_bmin,m_bmax,1.0e-5);
    k2=k1;
    k1=m_kfactor;
    m_kfactor=k2+(ratio-ratio1)*(m_kfactor-k2)/(ratio2-ratio1);
    msg_Debugging()<<"iterate r2 = "<<ratio2<<",\t r= "<<ratio<<",\t r2-r = "<<ratio2-ratio
      		   <<"\t => "<<m_kfactor<<"\t <- "<<k1<<std::endl;
    if (!(m_kfactor>0.0)) {
      msg_Error()<<"Profile_Function_Base::CalculateOMean("<<ratio<<"): "
			 <<"Cannot determine k."<<std::endl;
      delete gausso;
      delete gaussp;
      return false;
    }
  } while(ATOOLS::dabs(ratio2-ratio)>1.0e-4);
#ifdef DEBUG__Profile_Function_Base
  Test_Func_O testfunco(this);
  ATOOLS::Gauss_Integrator gaussh(&testfunco);
  double f_c=gaussh.Integrate(m_bmin,m_bmax,1.0e-5);
  double norm=gausso->Integrate(m_bmin,m_bmax,1.0e-5);
  double pmean=gaussp->Integrate(m_bmin,m_bmax,1.0e-5);
  double om=f_c;
  f_c/=norm;
  om/=pmean*m_kfactor;
  norm/=m_kfactor*om*pmean;
#endif
  delete gausso;
  delete gaussp;
  m_omean=ratio2/m_kfactor;
#ifdef DEBUG__Profile_Function_Base
  Test_Func_F testfuncf(this);
  ATOOLS::Gauss_Integrator gaussf(&testfuncf);
  double fmean=gaussf.Integrate(m_bmin,m_bmax,1.0e-5);
  fmean/=pmean;
#endif
  msg_Info()<<"Profile_Function_Base::CalculateOMean("<<ratio<<"): "
	    <<"Results are {\n   k           = "<<m_kfactor
#ifdef DEBUG__Profile_Function_Base
	    <<"\n   <\\tilde{O}> = "<<om<<" norm = "<<norm
	    <<"\n   f_c         = "<<f_c
	    <<"\n   <f(b)>      = "<<fmean<<" -> "<<fmean/f_c 
#endif
	    <<"\n   <\\tilde{O}> = "<<m_omean<<"\n}"<<std::endl;
  return true;
}

double Profile_Function_Base::DiceImpactParameter() const
{
  double b=0.0;
  double maxintegral=MajorIntegral(m_bmin);
  do {
    b=InverseMajorIntegral(ATOOLS::ran.Get()*maxintegral);
  } while (Value(b)<=ATOOLS::ran.Get()*MajorValue(b));
  return b;
}

