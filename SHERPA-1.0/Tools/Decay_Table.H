#ifndef Decay_Table_h
#define Decay_Table_h

#include "Flavour.H"

namespace SHERPA {
  class Decay_Channel {
    double               m_width, m_BR;
    ATOOLS::Flavour * fl;
  };
  /*!
    \class Decay_Channel
    \brief So far not being used

    This class is not used so far; work is ongoing.
  */
}

#endif
