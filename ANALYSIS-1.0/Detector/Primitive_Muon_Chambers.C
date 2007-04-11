#include "Primitive_Muon_Chambers.H"
#include "Random.H"
#include "MathTools.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Primitive_Muon_Chambers_Getter,"Muon_Chambers",Analysis_Object,Argument_Matrix);

Analysis_Object *				     
Primitive_Muon_Chambers_Getter::operator()(const Argument_Matrix & params) const  
{ 
  double   mineta(-5.), maxeta(5.), minphi(0.), maxphi(2.*M_PI);
  long int neta(1000000), nphi(628319);
  muonprob::code mp=muonprob::full;
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
    else if (cur[0]=="Muon_Mode" && cur.size()==2) {
      mp=muonprob::code(ToType<int>(cur[1]));
    }
    else if (cur[0]=="Miss_Probability" && cur.size()==2 && mp==muonprob::flatmiss) {
      missprob=ToType<double>(cur[1]);
    }
    else if (cur[0]=="Threshold" && cur.size()==2) {
      threshold=ToType<double>(cur[1]);
    }
    else if (cur[0]=="Resolution" && cur.size()==2) {
      resolution=ToType<double>(cur[1]);
    }
  }
  Primitive_Muon_Chambers * 
    muons(new Primitive_Muon_Chambers(neta,nphi,mineta,maxeta,minphi,maxphi));
  muons->SetMuonMode(mp);
  muons->SetThreshold(threshold);
  muons->SetResolution(resolution);
  muons->SetMissProbability(missprob);
  return muons; 
}

void Primitive_Muon_Chambers_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Muon_Chambers"; 
}

Primitive_Muon_Chambers::
Primitive_Muon_Chambers(const long int neta,const long int nphi,
			const double mineta,const double maxeta,
			const double minphi,const double maxphi) :
  Primitive_Detector_Element(pde::tracks,neta,nphi,mineta,maxeta,minphi,maxphi),
  m_muonprob(muonprob::full),
  m_threshold(0.), m_missprob(0.), m_resolution(0.)
{
  p_qualifier = Particle_Qualifier_Getter::GetObject("92","muon");
  m_name      = "MuonChambers";
}


Primitive_Muon_Chambers::~Primitive_Muon_Chambers() {}

Analysis_Object * Primitive_Muon_Chambers::GetCopy() const {
  Primitive_Muon_Chambers * muon_chambers =
    new Primitive_Muon_Chambers(m_neta,m_nphi,m_etamin,m_etamax,m_phimin,m_phimax);
  muon_chambers->SetMissProbability(m_missprob);
  muon_chambers->SetThreshold(m_threshold);
  muon_chambers->SetResolution(m_resolution);
  muon_chambers->SetMuonMode(m_muonprob);
  return muon_chambers;
}

void Primitive_Muon_Chambers::Fill(const Particle_List * plist) {
  Reset();
  for (Particle_List::const_iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!(*p_qualifier)((*pit))) continue;
    FillParticleInDetectorElement((*pit));
  }
}

bool Primitive_Muon_Chambers::FillParticleInDetectorElement(const Particle * part) {
  long int etapos,phipos;
  MatchCell(part->Momentum(),etapos,phipos);
  if (etapos<0||phipos<0) return false;
  std::pair<long int, long int> pos;
  pos.first=etapos; pos.second=phipos;
  double energydeposit = EnergyDeposit(part);
  if (m_tracks.find(pos)!=m_tracks.end()) 
    m_tracks.find(pos)->second += energydeposit;
  else m_tracks[pos] = energydeposit;
  return false;
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

void Primitive_Muon_Chambers::Print(std::ostream &) {}

void Primitive_Muon_Chambers::AddNoise() {}
