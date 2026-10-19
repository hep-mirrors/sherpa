#ifndef HADRONS_ME_Library_Exchange_ME_Shim_H
#define HADRONS_ME_Library_Exchange_ME_Shim_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "METOOLS/HadronCurrents/Exchange_ME.H"

namespace HADRONS {

  // The amplitude itself lives in METOOLS::Exchange_ME, so that EXTRA_XS and
  // the decay handler share ONE constructor and METOOLS stays the common
  // dependency of both. This is only the adaptor that makes it findable from
  // a decay channel.
  //
  // It exists because HD_ME_Base and METOOLS::ME_Base are the same class
  // written twice -- same members, same virtuals, both ultimately
  // METOOLS::Spin_Amplitudes -- but each DECLARE_GETTER declares its own
  // registry, so a class registered against one is invisible to lookups in
  // the other. A decay channel's ME block is resolved in the HD_ME_Base
  // registry.
  //
  // Collapsing the two bases would delete this file. That is a worthwhile
  // job -- 14 subclasses and 31 referencing files -- and deliberately not
  // this one, so that unifying the amplitude does not turn into migrating
  // HADRONS++.
  class Exchange_ME_Shim : public HD_ME_Base {
    METOOLS::Exchange_ME * p_me;
  public:
    Exchange_ME_Shim(const ATOOLS::Flavour_Vector& flavs,
                     const std::vector<int>& indices,
                     const std::string& name);
    ~Exchange_ME_Shim();
    void SetModelParameters(ATOOLS::Scoped_Settings& s);
    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };
}

#endif
