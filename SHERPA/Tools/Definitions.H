#ifndef SHERPA_Tools_Definitions_H
#define SHERPA_Tools_Definitions_H


namespace SHERPA {
  struct eventtype {
    enum code {
      StandardPerturbative = 1,
      MinimumBias          = 2,
      HadronDecay          = 3,
      NeutrinoNucleon      = 4,
      EventReader          = 10
    };
  };
}
#endif
