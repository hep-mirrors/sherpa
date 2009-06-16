#ifndef SHERPA_Tools_HepMC2_Interface_H
#define SHERPA_Tools_HepMC2_Interface_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"


#include <map>
#include <set>

namespace HepMC {
  class GenEvent;
  class GenVertex;
  class GenParticle;
}

namespace ATOOLS {
  class Blob_List;
  class Blob;
  class Particle;
}

namespace SHERPA {

  class HepMC2_Interface {
  private :

    std::map<ATOOLS::Blob*,HepMC::GenVertex *> m_blob2genvertex;
    std::map<ATOOLS::Particle*,HepMC::GenParticle *> m_particle2genparticle;

    HepMC::GenEvent *p_event;

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

    // inline functions
    inline HepMC::GenEvent *const GenEvent() { return p_event; }

  };// end of class HepMC2_Interface

}// end of namespace SHERPA

#endif
