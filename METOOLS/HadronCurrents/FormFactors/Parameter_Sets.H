#ifndef METOOLS_HadronCurrents_Parameter_Sets_H
#define METOOLS_HadronCurrents_Parameter_Sets_H

#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Phys/Flavour.H"
#include <vector>

namespace METOOLS {

  ///////////////////////////////////////////////////////////////////////////
  //
  // Local parameter sets.
  //
  // ParameterSets/*.yaml each record ONE published fit (provenance, quotes,
  // closure) plus a generated flat `sherpa:` block.  A set is selected by the
  // PAIR (channel, FORM_FACTOR) -- never by the model number alone, which many
  // channels share -- and is applied whole or not at all.
  //
  //   HADRON_DECAYS:
  //     Resonance_Parameters: global | local | local_else_global
  //     Parameter_Sets_Path:  <directory of the *.yaml files>
  //
  //   global             registry everywhere.  The default, and bit-identical
  //                      to running without this machinery.
  //   local              every form factor with a FORM_FACTOR must have a set
  //                      for its (channel, model); a missing one is fatal.
  //   local_else_global  a set where one exists, the registry otherwise -- per
  //                      CHANNEL, never per parameter.
  //
  // WHY THIS LIVES IN METOOLS.  A parameter set describes resonance poles and
  // widths, which is form-factor material; nothing about it is specific to a
  // hadron decay.  It was first written into HADRONS++'s Current_ME, because
  // that is the one place holding both the channel and the current, but the
  // cost was that only decays could use it -- e+e- production never passes
  // through Current_ME.  Applying it in FormFactor_Base instead reaches every
  // consumer of a form factor, and leaves HADRONS++ untouched.
  //
  ///////////////////////////////////////////////////////////////////////////

  //! Merge the set for (channel, ffmodel) into the model, per the policy above.
  //
  //  flavs is the current's full flavour vector and hadrons its decay indices.
  //  The channel is looked up first as every flavour but the first -- the
  //  convention the YAML files use, which for a decay is "the products",
  //  neutrino included -- and then as the hadrons alone, so that a production
  //  process whose first two entries are both beams can also be keyed.  Each
  //  is tried in its charge conjugate as well.
  //
  //  Does nothing when model is NULL or ffmodel names no dynamical model.
  void ApplyParameterSet(GeneralModel * model,
                         const ATOOLS::Flavour_Vector & flavs,
                         const std::vector<int> & hadrons,
                         const int & ffmodel);
}

#endif
