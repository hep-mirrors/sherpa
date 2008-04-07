#ifndef Soft_Photon_Handler_H
#define Soft_Photon_Handler_H

#include "Photons.H"
#include "Blob.H"
#include <string>

namespace SHERPA {
  struct softphotons {
    enum code {
      off        = 0,
      exp        = 1,
      exp_order1 = 2
    };
  };

  class Soft_Photon_Handler {
  private:
    softphotons::code   m_mode;
    PHOTONS::Photons  * p_yfs;
  public :
    // constructor
    Soft_Photon_Handler(std::string path=std::string("./"),
			std::string datfile=std::string("Fragmentation.dat"));
    // destructor
    ~Soft_Photon_Handler();
    bool AddRadiation(ATOOLS::Blob *);
  };// end of class Soft_Photon_Handler

}// end of namespace SHERPA

#endif

