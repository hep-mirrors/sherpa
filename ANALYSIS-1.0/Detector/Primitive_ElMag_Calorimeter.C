#include "Primitive_ElMag_Calorimeter.H"
#include "Random.H"
#include "MathTools.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Primitive_ElMag_Calorimeter_Getter,"ECal",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *				     
Primitive_ElMag_Calorimeter_Getter::operator()(const Argument_Matrix & params) const  
{ 
  double   mineta(-5.), maxeta(5.), minphi(0.), maxphi(2.*M_PI);
  long int neta(1000), nphi(628);
  emcalor::code calor=emcalor::full;
  emgeom::code geom=emgeom::full;
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
      calor=emcalor::code(ToType<int>(cur[1]));
    }
    else if (cur[0]=="Geometry_Mode" && cur.size()==2) {
      geom=emgeom::code(ToType<int>(cur[1]));
    }
    else if (cur[0]=="Miss_Probability" && cur.size()==2 && geom==emgeom::flatmiss) {
      missprob=ToType<double>(cur[1]);
    }
    else if (cur[0]=="Threshold" && cur.size()==2) {
      threshold=ToType<double>(cur[1]);
    }
    else if (cur[0]=="Resolution" && cur.size()==2) {
      resolution=ToType<double>(cur[1]);
    }
  }
  Primitive_ElMag_Calorimeter * 
    ecal(new Primitive_ElMag_Calorimeter(neta,nphi,mineta,maxeta,minphi,maxphi));
  ecal->SetCalorimetryMode(calor);
  ecal->SetGeometryMode(geom);
  ecal->SetThreshold(threshold);
  ecal->SetResolution(resolution);
  ecal->SetMissProbability(missprob);
  return ecal;
}

Primitive_ElMag_Calorimeter::
Primitive_ElMag_Calorimeter(const long int neta,const long int nphi,
			    const double mineta,const double maxeta,
			    const double minphi,const double maxphi) :
  Primitive_Detector_Element(pde::cells,neta,nphi,mineta,maxeta,minphi,maxphi),
  m_calor(emcalor::full),m_geom(emgeom::full),
  m_threshold(0.), m_missprob(0.), m_resolution(0.),
  p_electron(Particle_Qualifier_Getter::GetObject("91","electron")),
  p_photon(Particle_Qualifier_Getter::GetObject("22","photon")),
  p_muon(Particle_Qualifier_Getter::GetObject("92","muon")),
  p_charged_hadron(Particle_Qualifier_Getter::GetObject("1","charged hadron"))
{
  p_qualifier = NULL;
  m_name      = "ECal";
}


Primitive_ElMag_Calorimeter::~Primitive_ElMag_Calorimeter() {}

Analysis_Object * Primitive_ElMag_Calorimeter::GetCopy() const {
  Primitive_ElMag_Calorimeter * elmag_calorimeter =
    new Primitive_ElMag_Calorimeter(m_neta,m_nphi,m_etamin,m_etamax,m_phimin,m_phimax);
  elmag_calorimeter->SetMissProbability(m_missprob);
  elmag_calorimeter->SetThreshold(m_threshold);
  elmag_calorimeter->SetResolution(m_resolution);
  elmag_calorimeter->SetCalorimetryMode(m_calor);
  elmag_calorimeter->SetGeometryMode(m_geom);
  return elmag_calorimeter;
}

void Primitive_ElMag_Calorimeter::Reset()
{
  for (m_etastrip=m_cells.begin();m_etastrip!=m_cells.end();m_etastrip++) {
    for (std::vector<double>::iterator etaphi=(*m_etastrip).begin();
	 etaphi!=(*m_etastrip).end();etaphi++) 
      (*etaphi) = 0.;
  }
  AddNoise();
  for (std::vector<std::vector<int> >::iterator etas=m_id.begin();etas!=m_id.end();etas++) {
    for (std::vector<int>::iterator etaphi=(*etas).begin();etaphi!=(*etas).end();etaphi++)
      (*etaphi) = 0;
  }
}

int Primitive_ElMag_Calorimeter::GetFlav(const long int etapos, const long int phipos) const {
  return m_id[etapos][phipos];
}

void Primitive_ElMag_Calorimeter::Fill(const Particle_List * plist) {
  Reset();
  for (Particle_List::const_iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!(*p_electron)((*pit)) && !(*p_photon)((*pit)) &&
	!(*p_muon)((*pit)) && !(*p_charged_hadron)((*pit))) continue;
    FillParticleInDetectorElement((*pit));
  }
}

bool Primitive_ElMag_Calorimeter::FillParticleInDetectorElement(const Particle * part) {
  long int etapos,phipos;
  MatchCell(part->Momentum(),etapos,phipos);
  if (etapos<0||phipos<0) return false;
  bool eraseit(false);
  m_cells[etapos][phipos] += EnergyDeposit(part,eraseit);
  m_id[etapos][phipos] = part->Flav().IsAnti()?-int(part->Flav().Kfcode()):int(part->Flav().Kfcode());
  return eraseit;
}

double Primitive_ElMag_Calorimeter::EnergyDeposit(const Particle * part,bool & eraseit) {
  double E(part->Momentum()[0]);
  if (E<m_threshold) return 0.; 
  if ((*p_electron)(part)) {
    SmearEnergy(part->Flav(),E);
    GeometryEffects(part->Flav(),E,part->Momentum().Eta(),part->Momentum().Phi());
    eraseit = true;
  }
  return E;
}

void Primitive_ElMag_Calorimeter::SmearEnergy(const Flavour & flav,double & E) {
  switch (m_calor) {
    case emcalor::full:
    default:
      return;
  }
}

void Primitive_ElMag_Calorimeter::GeometryEffects(const Flavour & flav,double & E,
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

void Primitive_ElMag_Calorimeter::Print(std::ostream &) {}

void Primitive_ElMag_Calorimeter::AddNoise() {}
