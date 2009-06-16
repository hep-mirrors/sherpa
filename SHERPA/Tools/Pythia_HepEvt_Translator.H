#ifndef SHERPA_Tools_Pythia_HepEvt_Translator_H
#define SHERPA_Tools_Pythia_HepEvt_Translator_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include <map>


namespace SHERPA {
  class HepEvt_Interface;

  class Pythia_HepEvt_Translator {
  private:
    HepEvt_Interface                * p_interface;
    ATOOLS::Translation_Map           m_convertH2S;
    ATOOLS::Translation_Map::iterator m_piter, m_piter1;
    ATOOLS::Blob_List *               p_blobs;
    ATOOLS::Blob                    * p_signal;
    std::set<int>                     m_bunchints,m_signalints;
    std::set<int>::iterator           m_spiter;

    int      m_nhep, m_fl1, m_fl2, *p_isthep, *p_idhep, *p_jmohep, *p_jdahep;
    double   m_weight, m_Q2, m_x1, m_x2, *p_phep, *p_vhep;

    std::vector<int>            ueints;

    void CopyHepEvtBlock();
    void HepEvt2Particles();
    bool ReconstructBeamsAndBunches();
    bool ReconstructSignalBlob();
    bool ReconstructShowerBlob();
    bool ReconstructFragmentationBlob();
    void ReconstructDecayChain(ATOOLS::Blob *,const int);
    void CleanUp();
  public:
    Pythia_HepEvt_Translator(HepEvt_Interface *);
    ~Pythia_HepEvt_Translator() {}
    bool ConstructBlobs(ATOOLS::Blob_List * const);
  };
}

#endif
