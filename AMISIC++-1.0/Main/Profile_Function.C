#include "Profile_Function.H"
#include "MathTools.H"

using namespace AMISIC;

template <> Profile_Function_Base*
Profile_Function_Base::CreateProfile<Flat_Profile>(const std::string &type,
						   const std::vector<double> &parameters)
{
  if (type==std::string("Flat") && parameters.size()>0) {
    return new Flat_Profile(parameters[0]);
  }
  return NULL;
}

Flat_Profile::Flat_Profile(const double radius):
  Profile_Function_Base(pft::flat),
  m_radius(radius) {}

double Flat_Profile::Value(const double b) const
{
  return 1.0/ATOOLS::sqr(m_radius);
}

double Flat_Profile::Integral(const double b) const
{
  return M_PI;
}

template <> Profile_Function_Base*
Profile_Function_Base::CreateProfile<Exponential_Profile>(const std::string &type,
						   const std::vector<double> &parameters)
{
  if (type==std::string("Exponential") && parameters.size()>0) {
    return new Exponential_Profile(parameters[0]);
  }
  return NULL;
}

Exponential_Profile::Exponential_Profile(const double radius):
  Profile_Function_Base(pft::exponential),
  m_radius(radius) {}

double Exponential_Profile::Value(const double b) const
{
  return 0.5*exp(-b/m_radius)/ATOOLS::sqr(m_radius);
}

double Exponential_Profile::Integral(const double b) const
{
  return M_PI;
}

template <> Profile_Function_Base*
Profile_Function_Base::CreateProfile<Gaussian_Profile>(const std::string &type,
						   const std::vector<double> &parameters)
{
  if (type==std::string("Gaussian") && parameters.size()>0) {
    return new Gaussian_Profile(parameters[0]);
  }
  return NULL;
}

Gaussian_Profile::Gaussian_Profile(const double radius):
  Profile_Function_Base(pft::gaussian),
  m_radius(radius) {}

double Gaussian_Profile::Value(const double b) const
{
  return 0.5*exp(-0.5*ATOOLS::sqr(b/m_radius))/ATOOLS::sqr(m_radius);
}

double Gaussian_Profile::Integral(const double b) const
{
  return M_PI;
}

template <> Profile_Function_Base*
Profile_Function_Base::CreateProfile<Double_Gaussian_Profile>(const std::string &type,
						   const std::vector<double> &parameters)
{
  if (type==std::string("Double_Gaussian") && parameters.size()>2) {
    return new Double_Gaussian_Profile(parameters[0],parameters[1],parameters[2]);
  }
  return NULL;
}

Double_Gaussian_Profile::Double_Gaussian_Profile(const double radius1,
						 const double radius2,
						 const double partition):
  Profile_Function_Base(pft::double_gaussian),
  m_partition(partition)
{
  m_radius[0]=radius1;
  m_radius[1]=radius2;
}

double Double_Gaussian_Profile::Value(const double b) const
{
  double sumsqr=ATOOLS::sqr(m_radius[0])+ATOOLS::sqr(m_radius[1]);
  return 0.5*ATOOLS::sqr((1.-m_partition)/m_radius[0])*
    exp(-0.5*ATOOLS::sqr(b/m_radius[0]))+2.*m_partition*(1.-m_partition)/
    (sumsqr)*exp(-b*b/(sumsqr))+0.5*ATOOLS::sqr(m_partition/m_radius[1])*
    exp(-0.5*ATOOLS::sqr(b/m_radius[1]));
}

double Double_Gaussian_Profile::Integral(const double b) const
{
  return M_PI;
}

