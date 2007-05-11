#include "Reconstructed_Object.H"
#include "Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Reconstructed_Object::Reconstructed_Object(Flavour flav,const double E,
					   const double eta,const double phi) :
  m_flav(flav), m_E(E), m_eta(eta), m_phi(phi), m_mom(Vec4D(0.,0.,0.,0.)),
  p_cells(NULL), p_tracks(NULL) 
{ }

Reconstructed_Object::Reconstructed_Object(ATOOLS::Flavour flav,
					   std::vector<Cell *> * cells) :
  m_flav(flav), m_E(0.), m_eta(0.), m_phi(0.), m_mom(Vec4D(0.,0.,0.,0.)), 
  p_cells(cells), p_tracks(NULL) 
{ 
  for (std::vector<Cell *>::iterator cit=p_cells->begin();
       cit!=p_cells->end(); cit++) {
    m_mom += (*cit)->TotalDeposit()*(*cit)->Direction();
  }
  m_E   = m_mom[0];
  m_eta = m_mom.Eta();
  m_phi = m_mom.Phi();
}

Reconstructed_Object::~Reconstructed_Object() { 
  std::cout<<METHOD<<" for "<<this<<", cells = "<<p_cells<<", tracks = "<<p_tracks<<std::endl;
  if (p_cells)  { delete p_cells;  p_cells = NULL; }
  if (p_tracks) { delete p_tracks; p_tracks = NULL; }
}


Particle * Reconstructed_Object::CreateParticle() { 
  double costheta((exp(m_eta)-exp(-m_eta))/(exp(m_eta)+exp(-m_eta)));
  double sintheta(sqrt(1.-costheta*costheta));
  Vec3D  direction = Vec3D(cos(m_phi)*sintheta,sin(m_phi)*sintheta,costheta);
  double p2 = sqr(m_E)-sqr(m_flav.PSMass()), p=(p2<0)?m_E:sqrt(p2);
  Vec4D  mom = Vec4D(m_E,p*direction);
  return new Particle(0,m_flav,mom,'r');
}

bool Reconstructed_Object::IsIncluded(const Cell * cell) const {
  if (!p_cells || p_cells->size()==0) return false;
  for (size_t i=0;i<p_cells->size();i++) {
    if (cell==(*p_cells)[i]) return true;
  }
  return false;
}

bool Reconstructed_Object::IsIncluded(const Track * track) const {
  if (!p_tracks || p_tracks->size()==0) return false;
  for (size_t i=0;i<p_tracks->size();i++) {
    if (track==(*p_tracks)[i]) return true;
  }
  return false;
}
