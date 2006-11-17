#ifndef Herwig_HepEvt_Translator_h
#define Herwig_HepEvt_Translator_h

#include "Blob.H"
#include "Blob_List.H"
#include "Particle.H"
#include "Message.H"
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
