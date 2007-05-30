#include "QCD_Splitting_Kernels.H"

#include "MathTools.H"

using namespace PDF;
using namespace ATOOLS;

#define TR 0.5
#define CA 3.0
#define CF 4.0/3.0

Q_QG::Q_QG(const Flavour &q): 
  Splitting_Kernel(q,q,kf::gluon) {}

double Q_QG::operator()(const double &z) const
{
  return CF*(1.0+z*z)/(1.0-z);
}

double Q_QG::Integral(const double &zmin,const double &zmax) const
{
  return CF*(0.5*(sqr(1.0+zmin)-sqr(1.0+zmax))-2.0*log((1.0-zmax)/(1.0-zmin)));
}

Q_GQ::Q_GQ(const Flavour &q): 
  Splitting_Kernel(q,kf::gluon,q), m_swaped(q) {}

double Q_GQ::operator()(const double &z) const
{
  return m_swaped(1.0-z);
}

double Q_GQ::Integral(const double &zmin,const double &zmax) const
{
  return m_swaped.Integral(1.0-zmin,1.0-zmax);
}

G_QQ::G_QQ(const Flavour &q): 
  Splitting_Kernel(kf::gluon,q,q.Bar()) {}

double G_QQ::operator()(const double &z) const
{
  return TR*(sqr(z)+sqr(1.0-z));
}

double G_QQ::Integral(const double &zmin,const double &zmax) const
{
  return TR*(zmax*(zmax*(2.0*zmax/3.0-1.0)+1.0)+
	     -zmin*(zmin*(2.0*zmin/3.0-1.0)+1.0));
}

G_GG::G_GG(): 
  Splitting_Kernel(kf::gluon,kf::gluon,kf::gluon) {}

double G_GG::operator()(const double &z) const
{
  // due to symmetry this is not P_{gg}(z) but 0.5*P_{gg}(z)
  return CA*sqr(1.0-z*(1.0-z))/(z*(1.0-z));
}

double G_GG::Value(const double &z) const
{
  return 2.0*(*this)(z);
}

double G_GG::Integral(const double &zmin,const double &zmax) const
{
  return CA*(log(zmax*(1.-zmin)/(zmin*(1.-zmax)))-
	     zmax*(zmax*(zmax/3.-0.5)+2.)+zmin*(zmin*(zmin/3.-0.5)+2.));
}
