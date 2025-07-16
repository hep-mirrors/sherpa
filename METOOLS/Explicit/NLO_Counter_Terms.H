#ifndef METOOLS_Explicit_NLO_Counter_Terms_H
#define METOOLS_Explicit_NLO_Counter_Terms_H

#include <cstring>

namespace ATOOLS {
  class Flavour;
}

namespace PDF {
  class PDF_Base;
}

namespace MODEL {
  class One_Running_AlphaS;
}

namespace METOOLS {
  double Hab(const ATOOLS::Flavour& a, const ATOOLS::Flavour& b);
  double FPab(const ATOOLS::Flavour& a, const ATOOLS::Flavour& b,
              const double& z);
  double SPab(const ATOOLS::Flavour& a, const ATOOLS::Flavour& b,
              const double& z);
  double IPab(const ATOOLS::Flavour& a, const ATOOLS::Flavour& b,
              const double& x);

  double AlphaSCounterTerm(const double& muR2,const double& muR2ref,
                           const double& as,
                           MODEL::One_Running_AlphaS * oras,
                           const size_t oqcd);
  double CollinearCounterTerms(const ATOOLS::Flavour &fl,const double &x,
                               const double &z,const double& as,
                               const double &t1,const double &t2,
                               const double &tref,
                               PDF::PDF_Base * pdf);
}

#endif
