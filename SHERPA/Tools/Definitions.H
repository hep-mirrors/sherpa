#ifndef SHERPA_Tools_Definitions_H
#define SHERPA_Tools_Definitions_H


namespace SHERPA {
  struct eventtype {
    enum code {
      StandardPerturbative,
      MinimumBias,
      HadronDecay,
      EventReader
    };
  };
}
#endif
