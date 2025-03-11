#ifndef ALPACA_Main_Alpaca_H
#define ALPACA_Main_Alpaca_H

#include "ALPACA/Tools/Parton.H"
#include "ALPACA/Tools/P2P_Translator.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace ALPACA {

  class Alpaca {
  private:
    P2P_Translator m_translator;
    std::vector<Parton *> m_partons;

    bool Harvest(ATOOLS::Blob_List * blobs);
    bool AddBlob(ATOOLS::Blob_List * blobs);
    bool Evolve();
  public:
    //constructor
    Alpaca();

    // destructor
    ~Alpaca() {}

    bool operator()(ATOOLS::Blob_List * blobs);
  };
}

#endif
