#include "AddOns/Analysis/Detector/Muon_Chambers.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Muon_Chambers_Getter,"Muon_Chambers",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Muon_Chambers_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'

  Muon_Chambers * mc = new Muon_Chambers(parameters());
  double   etamin(0.),etamax(0.);

  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur.size()<2) continue;
    else if (cur[0]=="Acceptance") {
      etamin  = ATOOLS::ToType<double>(cur[1]);
      etamax  = ATOOLS::ToType<double>(cur[2]);
      mc->SetAcceptance(etamin,etamax);
    }
  }
  return mc;
}									

void Muon_Chambers_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Acceptance   Parameters=etamin, etamax"<<std::endl; 
}

Muon_Chambers::Muon_Chambers(Primitive_Analysis * ana) :
  Detector_Element(ana,"Muon_Chambers")
{
  m_isdet = true;
}

Muon_Chambers::~Muon_Chambers() { Reset(); }

void Muon_Chambers::SetAcceptance(const double etamin,const double etamax) {
  m_etamin = etamin; m_etamax = etamax;
}

Analysis_Object * Muon_Chambers::GetCopy() const {
  Muon_Chambers * mc = new Muon_Chambers(p_ana);
  mc->SetAcceptance(m_etamin,m_etamax);
  return mc;
}


void Muon_Chambers::Reset() {
  for (std::list<Track *>::iterator trit=m_tracks.begin();
       trit!=m_tracks.end();trit++) delete (*trit);
  m_tracks.clear();
}

bool Muon_Chambers::Fill(const double E,const double eta,
		   const double phi,ATOOLS::Particle * part) {
  if (eta>m_etamax || eta<m_etamin) return false;
  Track * track = new Track;
  track->eta  = eta;
  track->phi  = phi;
  track->mom  = part->Momentum();
  track->flav = part->Flav();
  m_tracks.push_back(track);
  return true;
}


