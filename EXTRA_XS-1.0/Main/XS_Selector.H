#ifndef XS_Selector_H
#define XS_Selector_H

#include "Single_XS.H"
#include "XS_Group.H"

namespace EXTRAXS { class XS_Group; }

namespace EXTRAXS {
  class XS_Selector {
  public:
    XS_Selector() {};
    static Single_XS * GetXS(int,int,ATOOLS::Flavour *);
    static bool FindInGroup(XS_Group *,XS_Base *&,int,int,ATOOLS::Flavour *);
  };
}

#endif
