#ifndef XS_Selector_H
#define XS_Selector_H

#include "Single_XS.H"
#include "XS_Group.H"

namespace EXTRAXS { class XS_Group; }

namespace EXTRAXS {
  class XS_Selector {
  public:
    XS_Selector() {};
    Single_XS * GetXS(int,int,APHYTOOLS::Flavour *);
    bool FindInGroup(XS_Group *,XS_Base *&,int,int,APHYTOOLS::Flavour *);
  };
}

#endif
