#ifndef HD_Selector_H
#define HD_Selector_H


#include "HD_ME_Base.H"
#include "Vector.H"
#include <string>

namespace HADRONS {
  class HD_ME_Selector {
    int m_mode;
    void SelectLightPseudoScalarDecay(int,ATOOLS::Flavour *,HD_ME_Base *&);
    void SetVector_For_2PS(int,int,int,ATOOLS::Flavour *,HD_ME_Base *);
    void SetDalitz_For_3P(int,ATOOLS::Flavour *,HD_ME_Base *);

    bool IsPseudoScalar(ATOOLS::Flavour &);
  public:
    HD_ME_Selector()      : m_mode(0) {};
    HD_ME_Selector(int m) : m_mode(m) {};
    HD_ME_Base * GetME(int,int,ATOOLS::Flavour *,
		       std::string _met=std::string(""));
  };

  inline bool HD_ME_Selector::IsPseudoScalar(ATOOLS::Flavour & flav) {
    int kfc = flav.Kfcode(); 
    if (kfc<1000 && kfc>110 && kfc%10==1) return true;
    return false;
  };

}

#endif
