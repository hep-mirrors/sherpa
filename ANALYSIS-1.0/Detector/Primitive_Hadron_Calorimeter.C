#include "Primitive_Hadron_Calorimeter.H"
#include "Random.H"
#include "MathTools.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Primitive_Hadron_Calorimeter_Getter,"HCal",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *				     
Primitive_Hadron_Calorimeter_Getter::operator()(const Argument_Matrix & params) const  
{ 
  double   mineta(-5.), maxeta(5.), minphi(0.), maxphi(2.*M_PI);
  long int neta(1000), nphi(628);
  hadcalor::code calor=hadcalor::full;
  hadgeom::code geom=hadgeom::full;
  double threshold(0.),missprob(0.),resolution(0.);
  for (size_t i=0;i<params.size();++i) {
    const std::vector<std::string> &cur=params[i];
    if (cur[0]=="Acceptance" && cur.size()==5) {
      mineta=ToType<double>(cur[1]);maxeta=ToType<double>(cur[2]);
      minphi=ToType<double>(cur[3]);maxphi=ToType<double>(cur[4]);
    }
    else if (cur[0]=="Granularity" && cur.size()==3) {
      neta=ToType<long int>(cur[1]);nphi=ToType<long int>(cur[2]);
    }
    else if (cur[0]=="Calorimetry_Mode" && cur.size()==2) {
      calor=hadcalor::code(ToType<int>(cur[1]));
    }
    else if (cur[0]=="Geometry_Mode" && cur.size()==2) {
      geom=hadgeom::code(ToType<int>(cur[1]));
    }
    else if (cur[0]=="Miss_Probability" && cur.size()==2 && geom==hadgeom::flatmiss) {
      missprob=ToType<double>(cur[1]);
    }
    else if (cur[0]=="Threshold" && cur.size()==2) {
      threshold=ToType<double>(cur[1]);
    }
    else if (cur[0]=="Resolution" && cur.size()==2) {
      resolution=ToType<double>(cur[1]);
    }
  }
  Primitive_Hadron_Calorimeter * 
    hcal(new Primitive_Hadron_Calorimeter(neta,nphi,mineta,maxeta,minphi,maxphi));
  hcal->SetCalorimetryMode(calor);
  hcal->SetGeometryMode(geom);
  hcal->SetThreshold(threshold);
  hcal->SetResolution(resolution);
  hcal->SetMissProbability(missprob);
  return hcal;
}

Primitive_Hadron_Calorimeter::
Primitive_Hadron_Calorimeter(const long int neta,const long int nphi,
			     const double mineta,const double maxeta,
			     const double minphi,const double maxphi) :
  Primitive_Detector_Element(pde::cells,neta,nphi,mineta,maxeta,minphi,maxphi),
  m_calor(hadcalor::full),m_geom(hadgeom::full),
  m_threshold(0.), m_missprob(0.), m_resolution(0.)
{
  p_qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject("hadron","Hadron"); 
  m_name      = "HCal";
}


Primitive_Hadron_Calorimeter::~Primitive_Hadron_Calorimeter() {}

Analysis_Object * Primitive_Hadron_Calorimeter::GetCopy() const {
  Primitive_Hadron_Calorimeter * hadron_calorimeter =
    new Primitive_Hadron_Calorimeter(m_neta,m_nphi,m_etamin,m_etamax,m_phimin,m_phimax);
  hadron_calorimeter->SetMissProbability(m_missprob);
  hadron_calorimeter->SetThreshold(m_threshold);
  hadron_calorimeter->SetResolution(m_resolution);
  hadron_calorimeter->SetCalorimetryMode(m_calor);
  hadron_calorimeter->SetGeometryMode(m_geom);
  return hadron_calorimeter;
}

void Primitive_Hadron_Calorimeter::Fill(const Particle_List * plist) {
  Reset();
  for (Particle_List::const_iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!(*p_qualifier)((*pit))) continue;
    FillParticleInDetectorElement((*pit));
  }
}

bool Primitive_Hadron_Calorimeter::FillParticleInDetectorElement(const Particle * part) {
  long int etapos,phipos;
  MatchCell(part->Momentum(),etapos,phipos);
  if (etapos<0||phipos<0) return false;
  m_cells[etapos][phipos] += EnergyDeposit(part);
  return true;
}

double Primitive_Hadron_Calorimeter::EnergyDeposit(const Particle * part) {
  double E(part->Momentum()[0]);
  if (E<m_threshold) return 0.; 
  SmearEnergy(part->Flav(),E);
  GeometryEffects(part->Flav(),E,part->Momentum().Eta(),part->Momentum().Phi());
  return E;
}

void Primitive_Hadron_Calorimeter::
SmearEnergy(const Flavour & flav,double & E) {
  switch (m_calor) {
    case hadcalor::full:
    default:
      return;
  }
}

void Primitive_Hadron_Calorimeter::
GeometryEffects(const Flavour & flav,double & E,const double eta,const double phi) {
  switch (m_geom) {
    case hadgeom::flatmiss:
      if (ran.Get()<m_missprob) { E=0.; return; }
    case hadgeom::full:
    default:
      return;
  }
}

void Primitive_Hadron_Calorimeter::Extract(Particle_List *) {}

void Primitive_Hadron_Calorimeter::Print(std::ostream &) {}

void Primitive_Hadron_Calorimeter::AddNoise() {}
