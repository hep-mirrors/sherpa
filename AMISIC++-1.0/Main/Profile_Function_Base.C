#include "Profile_Function_Base.H"

#include "Profile_Function.H"

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

Profile_Function_Base::Profile_Function_Base(const pft::code code):
  m_type(code) {}

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

