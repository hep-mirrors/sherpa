#include "Primitive_ElMag_Calorimeter.H"
#include "Random.H"
#include "MathTools.H"

using namespace ANALYSIS;
using namespace ATOOLS;


Primitive_ElMag_Calorimeter::
Primitive_ElMag_Calorimeter(const int neta,const int nphi,
			     const double mineta=-5.,const double maxeta=5.,
			     const double minphi=0.,const double maxphi=2.*M_PI) :
  Primitive_Detector_Element(neta,nphi,mineta,maxeta,minphi,maxphi),
  m_calor(emcalor::full),m_geom(emgeom::full),
  m_thres(0.), m_missprob(0.), m_resolution(0.),
  p_electron(Particle_Qualifier_Getter::GetObject("91","electron")),
  p_photon(Particle_Qualifier_Getter::GetObject("22","photon")),
  p_muon(Particle_Qualifier_Getter::GetObject("92","muon")),
  p_charged_hadron(Particle_Qualifier_Getter::GetObject("1","charged hadron"))
{
  p_qualifier = Particle_Qualifier_Getter::GetObject((("91","electron") | ("22","photon")) |
						     (("92","muon") | ("1","charged_hadron")));
  m_name      = "El-Mag Calorimeter";
}


Primitive_ElMag_Calorimeter::~Primitive_ElMag_Calorimeter() {}

void Primitive_ElMag_Calorimeter::Fill(const Particle_List * plist) {
  Reset();
  int etapos, phipos;
  for (Particle_List::iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!p_qualifier((*pit))) continue;
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
  if (p_electron(part)) {
  SmearEnergy(flav,E);
  GeometryEffects(flav,E,mom.Eta(),mom.Phi());
  return E;
}

void Primitive_ElMag_Calorimeter::
SmearEnergy(const Flavour & const flav,double & E) {
  switch (m_calor) {
    case emcalor::full:
    default:
      return;
  }
}

void Primitive_ElMag_Calorimeter::
GeometryEffects(const Flavour & const flav,double & E,
		const double eta,const double phi) {
  switch (m_geom) {
    case emgeom::flatmiss:
      if (ran.Get()<m_missprob) { E=0.; return; }
    case emgeom::full:
    default:
      return;
  }
}

void Primitive_ElMag_Calorimeter::Extract(Particle_List *);

void Primitive_ElMag_Calorimeter::Reset()
{
  for (int i=0; i<m_nx; ++i) {
    for (int j=0; j<m_ny; ++j) p_cells[i][j]=0.;
  }
  AddNoise();
}

void Print(std::ostream & =std::cout);

void Primitive_ElMag_Calorimeter::AddNoise() {}
