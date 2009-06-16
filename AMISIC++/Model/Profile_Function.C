#include "AMISIC++/Model/Profile_Function.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;

namespace AMISIC {

  template <> Profile_Function_Base*
  Profile_Function_Base::CreateProfile<Flat_Profile>(const std::string &type,
						     const std::vector<double> &parameters)
  {
    if (type==std::string("Flat") && parameters.size()>0) {
      return new Flat_Profile(parameters[0]);
    }
    return NULL;
  }

}

Flat_Profile::Flat_Profile(const double radius):
  Profile_Function_Base(pft::flat,0.0,radius),
  m_radius(radius) 
{
  m_omax=Value(m_bmin);
  m_omin=Value(m_bmax);
  m_norm=M_PI;
  THROW(not_implemented,"Flat profile not implemented yet");
}


double Flat_Profile::Value(const double b) const
{
  return 0.0;
}

double Flat_Profile::MajorValue(const double b) const
{
  return 0.0;
}

double Flat_Profile::MajorIntegral(const double b) const
{
  return 0.0;
}

double Flat_Profile::InverseMajorIntegral(const double I) const
{
  return 0.0;
}

namespace AMISIC {

  template <> Profile_Function_Base*
  Profile_Function_Base::CreateProfile<Exponential_Profile>(const std::string &type,
							    const std::vector<double> &parameters)
  {
    if (type==std::string("Exponential") && parameters.size()>0) {
      return new Exponential_Profile(parameters[0]);
    }
    return NULL;
  }

}

Exponential_Profile::Exponential_Profile(const double radius):
  Profile_Function_Base(pft::exponential,0.0,10.0*radius),
  m_radius(radius) 
{
  m_omax=Value(m_bmin);
  m_omin=Value(m_bmax);
  m_norm=M_PI;
}

double Exponential_Profile::Value(const double b) const
{
  return 0.5*exp(-b/m_radius)/m_radius;
}

double Exponential_Profile::MajorValue(const double b) const
{
  return 0.5*exp(-b/m_radius)/m_radius;
}

double Exponential_Profile::MajorIntegral(const double b) const
{
  return -M_PI*(exp(-b/m_radius)-exp(-m_bmax/m_radius));
}

double Exponential_Profile::InverseMajorIntegral(const double I) const
{
  return -m_radius*log(-I/M_PI+exp(-m_bmax/m_radius));
}

namespace AMISIC {

  template <> Profile_Function_Base*
  Profile_Function_Base::CreateProfile<Gaussian_Profile>(const std::string &type,
							 const std::vector<double> &parameters)
  {
    if (type==std::string("Gaussian") && parameters.size()>0) {
      return new Gaussian_Profile(parameters[0]);
    }
    return NULL;
  }

}

Gaussian_Profile::Gaussian_Profile(const double radius):
  Profile_Function_Base(pft::gaussian,0.0,10.0*radius),
  m_radius(radius) 
{
  m_omax=Value(m_bmin);
  m_omin=Value(m_bmax);
  m_norm=exp(-0.5*ATOOLS::sqr(m_bmin/m_radius));
  m_norm-=exp(-0.5*ATOOLS::sqr(m_bmax/m_radius));
  m_norm*=M_PI;
}

double Gaussian_Profile::Value(const double b) const
{
  return 0.5*exp(-0.5*ATOOLS::sqr(b/m_radius))/ATOOLS::sqr(m_radius);
}

double Gaussian_Profile::MajorValue(const double b) const
{
  return 0.5*exp(-0.5*ATOOLS::sqr(b/m_radius))/ATOOLS::sqr(m_radius);
}

double Gaussian_Profile::MajorIntegral(const double b) const
{
  return M_PI*(exp(-0.5*ATOOLS::sqr(b/m_radius))-
	       exp(-0.5*ATOOLS::sqr(m_bmax/m_radius)));
}

double Gaussian_Profile::InverseMajorIntegral(const double I) const
{
  return m_radius*sqrt(-2.0*log(I/M_PI+exp(-0.5*ATOOLS::sqr(m_bmax/m_radius))));
}

namespace AMISIC {

  template <> Profile_Function_Base*
  Profile_Function_Base::CreateProfile<Double_Gaussian_Profile>(const std::string &type,
								const std::vector<double> &parameters)
  {
    if (type==std::string("Double_Gaussian") && parameters.size()>2) {
      return new Double_Gaussian_Profile(parameters[0],parameters[1],parameters[2]);
    }
    return NULL;
  }

}

Double_Gaussian_Profile::Double_Gaussian_Profile(const double radius1,
						 const double radius2,
						 const double partition):
  Profile_Function_Base(pft::double_gaussian,0.0,10.0*radius1),
  m_partition(partition)
{
  m_radius[0]=radius1;
  m_radius[1]=radius2;
  m_omax=Value(m_bmin);
  m_omin=Value(m_bmax);
  m_norm=M_PI;
  m_rmin=ATOOLS::Min(m_radius[0],m_radius[1]);
  m_rmax=ATOOLS::Max(m_radius[0],m_radius[1]);
}

double Double_Gaussian_Profile::Value(const double b) const
{
  double sumsqr=ATOOLS::sqr(m_radius[0])+ATOOLS::sqr(m_radius[1]);
  return 0.5*ATOOLS::sqr((1.-m_partition)/m_radius[0])*
    exp(-0.5*ATOOLS::sqr(b/m_radius[0]))+2.*m_partition*(1.-m_partition)/
    (sumsqr)*exp(-b*b/(sumsqr))+0.5*ATOOLS::sqr(m_partition/m_radius[1])*
    exp(-0.5*ATOOLS::sqr(b/m_radius[1]));
}

double Double_Gaussian_Profile::MajorValue(const double b) const
{
  return 0.5*exp(-0.5*ATOOLS::sqr(b/m_rmax))/ATOOLS::sqr(m_rmin);
}

double Double_Gaussian_Profile::MajorIntegral(const double b) const
{
  return M_PI*ATOOLS::sqr(m_rmax/m_rmin)*
    (exp(-0.5*ATOOLS::sqr(b/m_rmax))-exp(-0.5*ATOOLS::sqr(m_bmax/m_rmax)));
}

double Double_Gaussian_Profile::InverseMajorIntegral(const double I) const
{
  return m_rmax*sqrt(-2.0*log(ATOOLS::sqr(m_rmin/m_rmax)*I/M_PI+
			      exp(-0.5*ATOOLS::sqr(m_bmax/m_rmax))));
}

