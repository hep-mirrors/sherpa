#ifndef Decay_Table_h
#define Decay_Table_h

#include "Flavour.H"

namespace SHERPA {
  class Decay_Channel {
    double               m_width, m_BR;
    APHYTOOLS::Flavour * fl;
  };
}

#endif
