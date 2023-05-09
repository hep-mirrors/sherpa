#ifndef NEUTRINOS_Tools_Propagator_Maps_H
#define NEUTRINOS_Tools_Propagator_Maps_H

#include "NEUTRINOS++/Tools/Propagator_Library.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <set>

namespace ATOOLS {
  class Scoped_Settings;
}

namespace NEUTRINOS {

  class Propagator_Maps {
    private:

      inline prop_type::code PropType(const std::string & form) {
      if (form=="none")             return prop_type::none;
      if (form=="massless")         return prop_type::massless;
      if (form=="massive")          return prop_type::massive;
      if (form=="unstable")         return prop_type::unstable;
      if (form=="unknown")          return prop_type::unknown;

      return prop_type::unknown;
    } 

    
    public:
      Propagator_Maps();
      ~Propagator_Maps();

      void Initialize();

      Propagator_Base * GetProp(kf_code & prop,prop_type::code & prop_type);
  };
  extern Propagator_Maps * ffprops;
}

#endif