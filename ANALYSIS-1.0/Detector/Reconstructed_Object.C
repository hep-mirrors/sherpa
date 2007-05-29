#include "Reconstructed_Object.H"
#include "Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Reconstructed_Object::Reconstructed_Object(Flavour flav,const double E,
					   const double eta,const double phi) :
  m_includetracks(false),
  m_flav(flav), m_E(E), m_eta(eta), m_phi(phi), m_mom(Vec4D(0.,0.,0.,0.))
{ 
}

Reconstructed_Object::Reconstructed_Object(Track * track) :
  m_includetracks(false),
  m_flav(track->flav), m_E(track->mom[0]), 
  m_eta(track->eta), m_phi(track->phi), m_mom(track->mom) 
{ 
  AddTrack(track);
}

Reconstructed_Object::Reconstructed_Object(ATOOLS::Flavour flav,
					   std::vector<Cell *> & cells) :
  m_includetracks(false),
  m_flav(flav), m_E(0.), m_eta(0.), m_phi(0.), m_mom(Vec4D(0.,0.,0.,0.))
{ 
  SetCells(cells);
  m_E   = m_mom[0];
  m_eta = m_mom.Eta();
  m_phi = m_mom.Phi();
}

Reconstructed_Object::~Reconstructed_Object() { 
  m_cells.clear();  
  m_tracks.clear(); 
}


void Reconstructed_Object::Update() {
  m_mom = Vec4D(0.,0.,0.,0.); 
  if (m_includetracks) {
    for (size_t i=0;i<m_tracks.size();i++) 
      m_mom += m_tracks[i]->mom;
  }
  for (size_t i=0;i<m_cells.size();i++)  
    m_mom += m_cells[i]->TotalDeposit()*m_cells[i]->Direction();

  m_E   = m_mom[0];
  m_eta = m_mom.Eta();
  m_phi = m_mom.Phi();
}

Vec4D Reconstructed_Object::TrueMom() const {
  Vec4D truemom(0.,0.,0.,0.),trackmom(0.,0.,0.,0.),cellmom(0.,0.,0.,0.); 
  if (m_includetracks) {
    for (size_t i=0;i<m_tracks.size();i++) trackmom += m_tracks[i]->mom;
  }
  std::set<Particle *> usedparts;
  std::map<ATOOLS::Particle *,double> * parts;
  Particle * part;
  for (size_t i=0;i<m_cells.size();i++) {
    parts = m_cells[i]->ParticleEntries();
    for (std::map<ATOOLS::Particle *,double>::iterator pit=parts->begin();
	 pit!=parts->end();pit++) {
      part = pit->first;
      if (usedparts.find(part->OriginalPart())==usedparts.end()) {
	cellmom += part->Momentum();
	usedparts.insert(part->OriginalPart());
      }
    }
  }
  return trackmom+cellmom;
}

Particle * Reconstructed_Object::CreateParticle() { 
  Update();
  Particle * part(new Particle(0,m_flav,m_mom,'r'));
  return part;
}

void Reconstructed_Object::SetUsed(const bool used) {
  for (size_t i=0;i<m_tracks.size();i++) m_tracks[i]->used = used;
  for (size_t i=0;i<m_cells.size();i++)  m_cells[i]->SetUsed(used);
}

bool Reconstructed_Object::IsIncluded(const Cell * cell) const {
  if (m_cells.size()==0) return false;
  for (size_t i=0;i<m_cells.size();i++) {
    if (cell==m_cells[i]) return true;
  }
  return false;
}

bool Reconstructed_Object::IsIncluded(const Track * track) const {
  if (m_tracks.size()==0) return false;
  for (size_t i=0;i<m_tracks.size();i++) {
    if (track==m_tracks[i]) return true;
  }
  return false;
}

void Reconstructed_Object::SetCells(const std::vector<Cell *> & cells) { 
  for (std::vector<Cell *>::const_iterator cit=cells.begin();
       cit!=cells.end(); cit++) {
    m_cells.push_back(*cit);
    m_mom += (*cit)->TotalDeposit()*(*cit)->Direction();
  }
}

void Reconstructed_Object::AddCell(Cell * cell) { m_cells.push_back(cell); }

void Reconstructed_Object::SetTracks(const std::vector<Track *> & tracks) { 
  for (std::vector<Track *>::const_iterator trit=tracks.begin();
       trit!=tracks.end(); trit++) {
    m_tracks.push_back(*trit);
    m_mom += (*trit)->mom;
  }
}

void Reconstructed_Object::AddTrack(Track * track) { m_tracks.push_back(track); }

int Reconstructed_Object::NumberOfParticles() const {
  int nop = m_tracks.size();
  for (size_t i=0;i<m_cells.size();i++) nop += m_cells[i]->ParticleEntries()->size();
  return nop;
}
