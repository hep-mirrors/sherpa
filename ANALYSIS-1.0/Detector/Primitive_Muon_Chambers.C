#include "Primitive_Muon_Chambers.H"
#include "Random.H"
#include "MathTools.H"

using namespace ANALYSIS;
using namespace ATOOLS;


Primitive_Muon_Chambers::
Primitive_Muon_Chambers(const int neta,const int nphi,
			     const double mineta,const double maxeta,
			     const double minphi,const double maxphi) :
  Primitive_Detector_Element(neta,nphi,mineta,maxeta,minphi,maxphi),
  m_muonprob(muonprob::full),
  m_threshold(0.), m_missprob(0.), m_resolution(0.)
{
  p_qualifier = Particle_Qualifier_Getter::GetObject("92","muon");
  m_name      = "Muon Chambers";
}


Primitive_Muon_Chambers::~Primitive_Muon_Chambers() {}

void Primitive_Muon_Chambers::Fill(const Particle_List * plist) {
  Reset();
  for (Particle_List::const_iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!(*p_qualifier)((*pit))) continue;
    FillParticleInDetectorElement((*pit));
  }
}

void Primitive_Muon_Chambers::
FillParticleInDetectorElement(const Particle * part) {
  int etapos,phipos;
  MatchCell(part->Momentum(),etapos,phipos);
  if (etapos<0||phipos<0) return;
  p_cells[etapos][phipos] += EnergyDeposit(part);
}

double Primitive_Muon_Chambers::
EnergyDeposit(const Particle * part) {
  double E(part->Momentum()[0]);
  if (E<m_threshold) return 0.; 
  SmearEnergy(part->Flav(),E);
  GeometryEffects(part->Flav(),E,part->Momentum().Eta(),part->Momentum().Phi());
  return E;
}

void Primitive_Muon_Chambers::
SmearEnergy(const Flavour & flav,double & E) {
}

void Primitive_Muon_Chambers::
GeometryEffects(const Flavour & flav,double & E,
		const double eta,const double phi) {
  switch (m_muonprob) {
    case muonprob::flatmiss:
      if (ran.Get()<m_missprob) { E=0.; return; }
    case muonprob::full:
    default:
      return;
  }
}

void Primitive_Muon_Chambers::Extract(Particle_List *) {}

void Primitive_Muon_Chambers::Reset()
{
  for (int i=0; i<m_neta; ++i) {
    for (int j=0; j<m_nphi; ++j) p_cells[i][j]=0.;
  }
  AddNoise();
}

void Primitive_Muon_Chambers::Print(std::ostream &) {}

void Primitive_Muon_Chambers::AddNoise() {}
