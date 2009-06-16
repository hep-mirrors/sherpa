#include "AddOns/Analysis/Detector/BTagging.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"

using namespace ANALYSIS;
using namespace ATOOLS;


BTagging::BTagging(const double tagprob,
		   const double mistaglight, const double mistagcharm) :
  m_tagprob(tagprob), m_mistaglight(mistaglight), m_mistagcharm(mistagcharm)
{}

BTagging::~BTagging() {}

bool BTagging::Tag(Reconstructed_Object * object) {
  if (FindB(object)) return TagIdentified();
  return Mistag();
}

bool BTagging::FindB(Reconstructed_Object * object) {
  m_foundcharm = false;
  std::vector<Cell *> cells(object->GetCells());
  if (!cells.empty()) {
    std::map<ATOOLS::Particle *,double> * parts(NULL);
    Particle * part(NULL);
    Blob     * prod(NULL);
    btp::code  btype;
    for (std::vector<Cell *>::iterator cit=cells.begin();
	 cit!=cells.end();cit++) {
      parts = (*cit)->ParticleEntries();
      if (!parts->empty()) {
	for (std::map<ATOOLS::Particle *,double>::iterator pit=parts->begin();
	     pit!=parts->end();pit++) {
	  part = pit->first->OriginalPart();
	  do {
	    prod  = part->ProductionBlob();
	    if (!prod) break;
	    btype = prod->Type(); 
	    part  = prod->InParticle(0);
	    if (part->Flav().IsB_Hadron()) return true;
	    if (part->Flav().IsC_Hadron()) m_foundcharm = true;
	  } while (btype!=btp::Fragmentation &&
		   btype!=btp::Cluster_Formation &&
		   btype!=btp::Cluster_Decay);
	}
      }
    }
  }
  return false;
}

bool BTagging::TagIdentified() {
  if (ran.Get()<m_tagprob) return true;
  return false;
}

bool BTagging::Mistag() {
  if (m_foundcharm && ran.Get()<m_mistagcharm) return true;
  if (ran.Get()<m_mistaglight) return true;
  return false;
}
