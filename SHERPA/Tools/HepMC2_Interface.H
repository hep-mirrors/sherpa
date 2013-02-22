#ifndef SHERPA_Tools_HepMC2_Interface_H
#define SHERPA_Tools_HepMC2_Interface_H
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__HEPMC2

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"


#include <map>
#include <set>

namespace HepMC {
  class GenEvent;
  class GenVertex;
  class GenParticle;
}

namespace ATOOLS {
  class Blob_List;
  class Particle;
}

namespace SHERPA {

  class HepMC2_Interface {
  private :

    std::map<ATOOLS::Blob*,HepMC::GenVertex *> m_blob2genvertex;
    std::map<ATOOLS::Particle*,HepMC::GenParticle *> m_particle2genparticle;
    std::set<ATOOLS::btp::code> m_ignoreblobs;

    HepMC::GenEvent *p_event;
    std::vector<HepMC::GenEvent*> m_subeventlist;

  public :

    // constructor
    HepMC2_Interface();

    // destructor
    ~HepMC2_Interface();

    // member functions
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs, HepMC::GenEvent& event,
                      double weight=1.0);
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs, double weight=1.0);
    bool Sherpa2HepMC(ATOOLS::Blob *,HepMC::GenVertex *&);
    bool Sherpa2HepMC(ATOOLS::Particle *,HepMC::GenParticle *&);

    bool Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                           HepMC::GenEvent& event,
                           double weight=1.0);
    bool Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs, double weight=1.0);

    // inline functions
    inline HepMC::GenEvent* GenEvent() { return p_event; }
    inline void Ignore(ATOOLS::btp::code type) { m_ignoreblobs.insert(type); }

    inline HepMC::GenEvent* GenSubEvent(size_t i)
    { if (i<m_subeventlist.size()) return m_subeventlist[i]; return NULL; }
    inline std::vector<HepMC::GenEvent*> GenSubEventList()
    { return m_subeventlist; }

    void DeleteGenSubEventList();

  };// end of class HepMC2_Interface

}// end of namespace SHERPA

#endif
#endif
