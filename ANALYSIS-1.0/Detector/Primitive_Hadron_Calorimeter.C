#include "Primitive_Hadron_Calorimeter.H"
#include "Random.H"
#include "MathTools.H"

using namespace ANALYSIS;
using namespace ATOOLS;


Primitive_Hadron_Calorimeter::
Primitive_Hadron_Calorimeter(const int neta,const int nphi,
			     const double mineta=-5.,const double maxeta=5.,
			     const double minphi=0.,const double maxphi=2.*M_PI) :
  Primitive_Detector_Element(neta,nphi,mineta,maxeta,minphi,maxphi),
  m_calor(hadcalor::full),m_geom(hadgeom::full),
  m_thres(0.), m_missprob(0.), m_resolution(0.)
{
  p_qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject("hadron","Hadron"); 
  m_name      = "Hadron Calorimeter";
}


Primitive_Hadron_Calorimeter::~Primitive_Hadron_Calorimeter() {}

void Primitive_Hadron_Calorimeter::Fill(const Particle_List * plist) {
  Reset();
  int etapos, phipos;
  for (Particle_List::iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!p_qualifier((*pit))) continue;
    FillParticleInDetectorElement((*pit));
  }
}

void Primitive_Hadron_Calorimeter::
FillParticleInDetectorElement(const Particle * part) {
  int etapos,phipos;
  MatchCell(part->Momentum(),etapos,phipos);
  if (etapos<0||phipos<0) return;
  p_cells[etapos][phipos] += EnergyDeposit(part->Flav(),part->Momentum());
}

double Primitive_Hadron_Calorimeter::
EnergyDeposit(const Particle * part) {
  double E(part->Momentum[0]);
  if (E<m_thres) return 0.; 
  SmearEnergy(part->Flav(),E);
  GeometryEffects(part->Flav(),E,mom.Eta(),mom.Phi());
  return E;
}

void Primitive_Hadron_Calorimeter::
SmearEnergy(const Flavour & const flav,double & E) {
  switch (m_calor) {
    case hadcalor::full:
    default:
      return;
  }
}

void Primitive_Hadron_Calorimeter::
GeometryEffects(const Flavour & const flav,double & E,const double eta,const double phi) {
  switch (m_geom) {
    case hadgeom::flatmiss:
      if (ran.Get()<m_missprob) { E=0.; return; }
    case hadgeom::full:
    default:
      return;
  }
}

void Primitive_Hadron_Calorimeter::Extract(Particle_List *) {}

void Primitive_Hadron_Calorimeter::Reset()
{
  for (int i=0; i<m_nx; ++i) {
    for (int j=0; j<m_ny; ++j) p_cells[i][j]=0.;
  }
  AddNoise();
}

void Print(std::ostream & =std::cout) {}

void Primitive_Hadron_Calorimeter::AddNoise() {}
