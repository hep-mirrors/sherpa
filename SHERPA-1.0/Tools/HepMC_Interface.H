#ifndef HepMC_Interface_h
#define HepMC_Interface_h

#include "Blob_List.H"
#include "Blob.H"
#include "Particle.H"
#include "Flavour.H"

#ifdef _USE_HEPMC_
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/ParticleDataTable.h"
#include "HepMC/ParticleData.h"
#else
namespace HepMC {
  class GenEvent;
  class ParticleDataTable;
  class ParticleData;
  class GenVertex;
  class GenParticle;
}
#endif

#include <map>
#include <set>

namespace SHERPA {
  typedef std::map<ATOOLS::Flavour,HepMC::ParticleData *>              Flavour2Particle;
  typedef std::map<ATOOLS::Flavour,HepMC::ParticleData *>::iterator    FlavsParts;
  typedef std::map<long int,HepMC::GenVertex *>                        Blob2Vertex;
  typedef std::map<long int,HepMC::GenVertex *>::iterator              BlobsVertices;
  typedef std::map<long int,HepMC::GenParticle *>                      Particle2Particle;
  typedef std::map<long int,HepMC::GenParticle *>::iterator            ParticlesParticles;

  class HepMC_Interface {
  private :
    HepMC::ParticleDataTable * p_particledatatable;
    Flavour2Particle           m_flavour2particle;
    Blob2Vertex                m_blob2vertex;
    Particle2Particle          m_particle2particle;

    void InitTheMap();
    bool Sherpa2HepMC(ATOOLS::Blob *,HepMC::GenVertex *&);
    bool Sherpa2HepMC(ATOOLS::Particle *,HepMC::GenParticle *&);
  public :
    HepMC_Interface();
    bool Sherpa2HepMC(ATOOLS::Blob_List *,HepMC::GenEvent *&); 
  };
}
#endif
