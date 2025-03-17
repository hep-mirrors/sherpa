#ifndef ALPACA_Tools_P2P_Translator_H
#define ALPACA_Tools_P2P_Translator_H

#include "ATOOLS/Phys/Particle.H"
#include "ALPACA/Tools/Parton.H"

namespace ALPACA {

  class P2P_Translator {
  private:
  public:
    //constructor
    P2P_Translator();

    // destructor
    ~P2P_Translator() {}

    Parton * operator()(ATOOLS::Particle * particle);
    ATOOLS::Particle * operator()(Parton * parton);
  };
}

#endif
