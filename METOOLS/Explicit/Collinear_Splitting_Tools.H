#ifndef METOOLS_Explicit_Collinear_Splitting_Tools_H
#define METOOLS_Explicit_Collinear_Splitting_Tools_H

namespace ATOOLS {
  class Flavour;
}

namespace METOOLS {
  double Hab(const ATOOLS::Flavour& a, const ATOOLS::Flavour& b);
  double FPab(const ATOOLS::Flavour& a, const ATOOLS::Flavour& b,
              const double& z);
  double SPab(const ATOOLS::Flavour& a, const ATOOLS::Flavour& b,
              const double& z);
  double IPab(const ATOOLS::Flavour& a, const ATOOLS::Flavour& b,
              const double& x);
}

#endif
