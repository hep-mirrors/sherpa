#ifndef HepMC_Interface_h
#define HepMC_Interface_h

#include "Blob_List.H"
#include "Blob.H"
#include "Particle.H"
#include "Flavour.H"

#ifdef CLHEP_SUPPORT
#include "CLHEP/HepMC/GenEvent.h"
#include "CLHEP/HepMC/GenVertex.h"
#include "CLHEP/HepMC/GenParticle.h"
#include "CLHEP/HepPDT/ParticleDataTable.hh"
#include "CLHEP/HepPDT/ParticleData.hh"

#include "CLHEP/HepPDT/DefaultConfig.hh"
#include "CLHEP/HepPDT/TableBuilder.hh"
#include "CLHEP/HepPDT/ParticleDataTableT.hh"
#include "CLHEP/HepPDT/TempParticleData.hh"

#else
namespace HepMC {
  class GenEvent;
  class GenVertex;
  class GenParticle;
}
namespace DefaultConfig {
  class ParticleDataTable;
}
namespace HepPDT {
  class ParticleData;
}
#endif

#include <map>
#include <set>

namespace SHERPA {
  typedef std::map<ATOOLS::Flavour,HepPDT::ParticleData *>           Flavour2Particle;
  typedef std::map<ATOOLS::Flavour,HepPDT::ParticleData *>::iterator FlavsParts;
  typedef std::map<long int,HepMC::GenVertex *>                        Blob2Vertex;
  typedef std::map<long int,HepMC::GenVertex *>::iterator              BlobsVertices;
  typedef std::map<long int,HepMC::GenParticle *>                      Parton2Particle;
  typedef std::map<long int,HepMC::GenParticle *>::iterator            PartonsParticles;

  class HepMC_Interface {
  private :
    DefaultConfig::ParticleDataTable * p_particledatatable;
    Blob2Vertex                m_blob2vertex;
    Parton2Particle            m_parton2particle;

    void InitTheMap();
    bool Sherpa2HepMC(ATOOLS::Blob *,HepMC::GenVertex *&);
    bool Sherpa2HepMC(ATOOLS::Particle *,HepMC::GenParticle *&);
  public :
    HepMC_Interface();
    ~HepMC_Interface();
    bool Sherpa2HepMC(ATOOLS::Blob_List *,HepMC::GenEvent *&); 
  };
}
#endif
