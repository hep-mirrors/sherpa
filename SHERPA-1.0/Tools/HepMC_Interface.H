#ifndef HepMC_Interface_H
#define HepMC_Interface_H

#include "CXXFLAGS.H"
#include "Blob_List.H"
#include "Particle.H"
#include "Message.H"

#include "CLHEP/HepMC/GenEvent.h"
#include "CLHEP/HepPDT/ParticleDataTable.hh"
#include "CLHEP/HepPDT/ParticleData.hh"
#include "CLHEP/HepPDT/DefaultConfig.hh"
#include "CLHEP/HepPDT/TableBuilder.hh"
#include "CLHEP/HepPDT/ParticleDataTableT.hh"
#include "CLHEP/HepPDT/TempParticleData.hh"

// namespace HepMC {
//   class GenEvent;
//   class GenVertex;
//   class GenParticle;
// }
// namespace DefaultConfig { class ParticleDataTable; }
// namespace HepPDT        { class ParticleData;      }


#include <map>
#include <set>

namespace SHERPA {

  typedef std::map<ATOOLS::Flavour,HepPDT::ParticleData *> Flavour2Particle;
  typedef std::map<long int,HepMC::GenVertex *>            Blob2Vertex;
  typedef std::map<long int,HepMC::GenParticle *>          Parton2Particle;

  class HepMC_Interface {
  private :

    DefaultConfig::ParticleDataTable * p_particledatatable;
    Blob2Vertex      m_blob2vertex;
    Parton2Particle  m_parton2particle;

    HepMC::GenEvent *p_event;
    bool             m_converted;

    void InitTheMap();
    bool Sherpa2HepMC(ATOOLS::Blob *,HepMC::GenVertex *&);
    bool Sherpa2HepMC(ATOOLS::Particle *,HepMC::GenParticle *&);

  public :

    // constructor
    HepMC_Interface();

    // destructor
    ~HepMC_Interface();

    // member functions
    void PrintEvent(int mode=1, std::ostream& ostr=ATOOLS::msg.Out());
    bool Sherpa2HepMC(ATOOLS::Blob_List *const blobs);
    inline void Reset() { m_converted = false; }

    // inline functions
    inline HepMC::GenEvent *const GenEvent() { return p_event; }

  };// end of class HepMC_Interface

}// end of namespace SHERPA

#endif
