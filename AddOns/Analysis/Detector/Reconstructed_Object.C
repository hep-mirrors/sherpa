#include "AddOns/Analysis/Detector/Reconstructed_Object.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Reconstructed_Object::Reconstructed_Object(Flavour flav,const double E,
					   const double eta,const double phi) :
  m_includetracks(false),
  m_flav(flav), m_E(E), m_eta(eta), m_phi(phi), 
  m_E_correction(0.), m_ET_correction(0.),
  m_mom(Vec4D(0.,0.,0.,0.))
{ 
}

Reconstructed_Object::Reconstructed_Object(Track * track) :
  m_includetracks(false),
  m_flav(track->flav), m_E(track->mom[0]), m_eta(track->eta), m_phi(track->phi), 
  m_E_correction(0.), m_ET_correction(0.),
  m_mom(track->mom) 
{ 
  AddTrack(track);
}

Reconstructed_Object::Reconstructed_Object(ATOOLS::Flavour flav,
					   std::vector<Cell *> & cells) :
  m_includetracks(false),
  m_flav(flav), m_E(0.), m_eta(0.), m_phi(0.), 
  m_E_correction(0.), m_ET_correction(0.),
  m_mom(Vec4D(0.,0.,0.,0.))
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
  Vec4D truemom(0.,0.,0.,0.);
  std::set<Particle *> parts;
  Particle * part;

  for (size_t i=0;i<m_cells.size();i++) {
    for (std::map<ATOOLS::Particle *,double>::iterator 
	   pit=m_cells[i]->ParticleEntries()->begin();
	 pit!=m_cells[i]->ParticleEntries()->end();pit++) {
      part = pit->first->OriginalPart();
      if (parts.find(part)==parts.end()) {
	truemom += part->Momentum();
	parts.insert(part);
      }
    }
  }
  if (m_includetracks) {
    for (size_t i=0;i<m_tracks.size();i++) {
      truemom += m_tracks[i]->mom;
    }
  }
  return truemom;
}

void Reconstructed_Object::CorrectTruth(const double val) {
  Vec4D truemom(0.,0.,0.,0.),depmom(0.,0.,0.,0.);
  std::set<Particle *> parts;
  Particle * part;

  for (size_t i=0;i<m_cells.size();i++) {
    depmom  += m_cells[i]->TotalDeposit()*m_cells[i]->Direction();
    for (std::map<ATOOLS::Particle *,double>::iterator 
	   pit=m_cells[i]->ParticleEntries()->begin();
	 pit!=m_cells[i]->ParticleEntries()->end();pit++) {
      part = pit->first->OriginalPart();
      if (parts.find(part)==parts.end()) {
	truemom += part->Momentum();
	parts.insert(part);
      }
    }
  }
  msg_Debugging()<<"   E correction : deposed = "<<depmom<<" vs. true = "<<truemom<<"."<<std::endl;
  double scaleit(truemom[0]/depmom[0]);
  double rana,dummy,factor;
  if (m_includetracks) {
    for (size_t i=0;i<m_tracks.size();i++) {
      do { ran.Gaussian(rana,dummy); } while (dabs(rana)>2.*M_PI);
      factor = 1.+val*rana/M_PI;
      m_tracks[i]->mom = factor*m_tracks[i]->mom;
    }
  }
  for (size_t i=0;i<m_cells.size();i++)  {
    do { ran.Gaussian(rana,dummy); } while (dabs(rana)>2.*M_PI);
    factor = scaleit*(1.+val*rana/M_PI);
    m_cells[i]->MultiplyDeposit(factor);
  }

  m_E_correction  = scaleit;
  m_ET_correction = truemom.EPerp()/depmom.EPerp(); 
}

void Reconstructed_Object::CorrectE(const double val) {
  if (m_includetracks) {
    for (size_t i=0;i<m_tracks.size();i++) m_tracks[i]->mom = val*m_tracks[i]->mom;
  }
  for (size_t i=0;i<m_cells.size();i++)  m_cells[i]->MultiplyDeposit(val);
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
