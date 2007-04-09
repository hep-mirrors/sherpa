#include "Primitive_ElMag_Calorimeter.H"
#include "Random.H"
#include "MathTools.H"

using namespace ANALYSIS;
using namespace ATOOLS;


Primitive_ElMag_Calorimeter::
Primitive_ElMag_Calorimeter(const int neta,const int nphi,
			     const double mineta,const double maxeta,
			     const double minphi,const double maxphi) :
  Primitive_Detector_Element(neta,nphi,mineta,maxeta,minphi,maxphi),
  m_calor(emcalor::full),m_geom(emgeom::full),
  m_thres(0.), m_missprob(0.), m_resolution(0.),
  p_electron(Particle_Qualifier_Getter::GetObject("91","electron")),
  p_photon(Particle_Qualifier_Getter::GetObject("22","photon")),
  p_muon(Particle_Qualifier_Getter::GetObject("92","muon")),
  p_charged_hadron(Particle_Qualifier_Getter::GetObject("1","charged hadron"))
{
  p_qualifier = NULL;
  m_name      = "El-Mag Calorimeter";
}


Primitive_ElMag_Calorimeter::~Primitive_ElMag_Calorimeter() {}

void Primitive_ElMag_Calorimeter::Fill(const Particle_List * plist) {
  Reset();
  for (Particle_List::const_iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!(*p_electron)((*pit)) && !(*p_photon)((*pit)) &&
	!(*p_muon)((*pit)) && !(*p_charged_hadron)((*pit))) continue;
    FillParticleInDetectorElement((*pit));
  }
}

void Primitive_ElMag_Calorimeter::
FillParticleInDetectorElement(const Particle * part) {
  int etapos,phipos;
  MatchCell(part->Momentum(),etapos,phipos);
  if (etapos<0||phipos<0) return;
  p_cells[etapos][phipos] += EnergyDeposit(part);
}

double Primitive_ElMag_Calorimeter::
EnergyDeposit(const Particle * part) {
  double E(part->Momentum()[0]);
  if (E<m_thres) return 0.; 
  if ((*p_electron)(part)) {
    SmearEnergy(part->Flav(),E);
    GeometryEffects(part->Flav(),E,part->Momentum().Eta(),part->Momentum().Phi());
  }
  return E;
}

void Primitive_ElMag_Calorimeter::
SmearEnergy(const Flavour & flav,double & E) {
  switch (m_calor) {
    case emcalor::full:
    default:
      return;
  }
}

void Primitive_ElMag_Calorimeter::
GeometryEffects(const Flavour & flav,double & E,
		const double eta,const double phi) {
  switch (m_geom) {
    case emgeom::flatmiss:
      if (ran.Get()<m_missprob) { E=0.; return; }
    case emgeom::full:
    default:
      return;
  }
}

void Primitive_ElMag_Calorimeter::Extract(Particle_List *) {}

void Primitive_ElMag_Calorimeter::Reset()
{
  for (int i=0; i<m_neta; ++i) {
    for (int j=0; j<m_nphi; ++j) p_cells[i][j]=0.;
  }
  AddNoise();
}

void Primitive_ElMag_Calorimeter::Print(std::ostream &) {}

void Primitive_ElMag_Calorimeter::AddNoise() {}
