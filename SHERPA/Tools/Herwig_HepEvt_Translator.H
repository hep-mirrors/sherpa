#ifndef SHERPA_Tools_Herwig_HepEvt_Translator_H
#define SHERPA_Tools_Herwig_HepEvt_Translator_H

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include <map>

namespace SHERPA {
  class HepEvt_Interface;

  class Herwig_HepEvt_Translator {
  private:
    HepEvt_Interface * p_interface;

    ATOOLS::Blob     * blob, * signal, * productionblob, * decayblob;
    ATOOLS::Particle * part, * mother, * daughter, * partner;

    std::vector<ATOOLS::Blob *> _blobs;
    std::set<int>               signalints, intermeds, outs, ins;
    std::set<int>::iterator     sit;
    std::vector<int>            ueints;
    ATOOLS::Translation_Map::iterator   piter, miter, partiter, diter;


    void FollowDaughters(const int &,ATOOLS::Blob *,ATOOLS::Blob *);
    void DeleteObsolete(const int);
    bool IdentifyBlobs(ATOOLS::Blob_List * const); 
    void HepEvt2Particle(const int);

  public:
    Herwig_HepEvt_Translator(HepEvt_Interface *);
    bool ConstructBlobs(ATOOLS::Blob_List * const);
  };
}

#endif
