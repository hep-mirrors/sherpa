#ifndef ALPACA_Main_Alpaca_H
#define ALPACA_Main_Alpaca_H

#include "ALPACA/EventGeneration/Event_Generator.H"
#include "ALPACA/Tools/Parton.H"
#include "ALPACA/Tools/P2P_Translator.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ALPACA/Tools/HeavyIon_Parameters.H"


namespace ALPACA {

  class Alpaca {
  private:
    P2P_Translator                                        m_translator;
    std::vector<ATOOLS::Particle *>                       m_inparticles;
    std::shared_ptr<std::list<std::shared_ptr<Parton>>>   p_partons;
    std::shared_ptr<Event_Generator>                      p_generator;
    std::shared_ptr<std::list<std::pair<int,int>>>        p_flow_backtrack;
    //MCEKRT_Handler                                      m_mcekrt; //MC-EKRT, temporarily commented out

    bool Harvest(ATOOLS::Blob_List * blobs);
    bool AddBlob(ATOOLS::Blob_List * blobs);
    bool Evolve();

    void CheckIfSinglet();
    void CheckFlowBacktrack(ATOOLS::Blob_List * blobs);
  public:
    //constructor
    Alpaca();

    // destructor
    ~Alpaca() {}

    bool operator()(ATOOLS::Blob_List * blobs);
  };
}

#endif
