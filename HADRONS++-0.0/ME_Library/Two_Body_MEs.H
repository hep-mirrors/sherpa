#ifndef Two_Body_MEs_H
#define Two_Body_MEs_H

#include "HD_ME_Base.H"

namespace HADRONS {
  class P_2Gamma : public HD_ME_Base {
  public:
    P_2Gamma(int,ATOOLS::Flavour *);
    double operator()(const ATOOLS::Vec4D *);
  };
}


#endif
