#ifndef Process_XS_Trait_H
#define Process_XS_Trait_H

#include "Process_Base.H"
#include "XS_Base.H"

namespace EXTRAXS {

  class Process_XS_Trait : AMEGIC::Process_Base {
    XS_Base * xs;
  public:
    Process_XS_Trait(XS_Base *);
  };

}



#endif
