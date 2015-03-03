#include "METOOLS/Explicit/Collinear_Splitting_Tools.H"

#include "ATOOLS/Phys/Flavour.H"

using namespace METOOLS;
using namespace ATOOLS;

double METOOLS::Hab(const Flavour &a,const Flavour &b)
{
  if (a.IsQuark()) {
    if (b.IsQuark()) return a==b?4.0/3.0*3.0/2.0:0.0;
    return 0.0;
  }
  else {
    if (b.IsQuark()) return 0.0;
    return 11.0/6.0*3.0-2.0/3.0*0.5*(Flavour(kf_jet).Size()/2);
  }
}

double METOOLS::FPab(const Flavour &a,const Flavour &b,const double &z)
{
  if (a.IsQuark()) {
    if (b.IsQuark()) return a==b?-4.0/3.0*(1.0+z):0.0;
    return 4.0/3.0*(1.0+sqr(1.0-z))/z;
  }
  else {
    if (b.IsQuark()) return 1.0/2.0*(z*z+sqr(1.0-z));
    return 3.0*2.0*((1.0-z)/z-1.0+z*(1.0-z));
  }
}

double METOOLS::SPab(const Flavour &a,const Flavour &b,const double &z)
{
  if (a.IsQuark()) {
    if (b.IsQuark()) return a==b?4.0/3.0*2.0/(1.0-z):0.0;
    return 0.0;
  }
  else {
    if (b.IsQuark()) return 0.0;
    return 3.0*2.0/(1.0-z);
  }
}

double METOOLS::IPab(const Flavour &a,const Flavour &b,const double &x)
{
  if (a.IsQuark()) {
    if (b.IsQuark() && a==b)
      return 4.0/3.0*2.0*log(1.0/(1.0-x));
    return 0.0;
  }
  else {
    if (b.IsQuark()) return 0.0;
    return 3.0*2.0*log(1.0/(1.0-x));
  }
}

