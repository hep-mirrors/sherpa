#ifndef SHERPA_PerturbativePhysics_Hard_Decay_Handler_H
#define SHERPA_PerturbativePhysics_Hard_Decay_Handler_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {
  class Decay_Map;
  class Decay_Table;
  class Decay_Channel;
  class Particle;
  class Blob;
}

namespace SHERPA {

  class Hard_Decay_Handler {
  private :
    ATOOLS::Decay_Map *      p_decaymap;
  public :
    Hard_Decay_Handler(std::string,std::string);
    ~Hard_Decay_Handler();

    void InitializeDecayMap();
  };
}

#endif

