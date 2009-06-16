#include "AddOns/Analysis/Detector/Tracker.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Tracker_Getter,"Tracker",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Tracker_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'

  Tracker * tracker = new Tracker(parameters());
  double   etamin(0.),etamax(0.);

  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur.size()<2) continue;
    else if (cur[0]=="Acceptance") {
      etamin  = ATOOLS::ToType<double>(cur[1]);
      etamax  = ATOOLS::ToType<double>(cur[2]);
      tracker->SetAcceptance(etamin,etamax);
    }
  }
  return tracker;
}									

void Tracker_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Acceptance   Parameters=etamin, etamax"<<std::endl; 
}

Tracker::Tracker(Primitive_Analysis * ana) :
  Detector_Element(ana,"Tracker")
{
  m_isdet = true;
}

Tracker::~Tracker() { Reset(); }

void Tracker::SetAcceptance(const double etamin,const double etamax) {
  m_etamin = etamin; m_etamax = etamax;
}

Analysis_Object * Tracker::GetCopy() const {
  Tracker * tracker = new Tracker(p_ana);
  tracker->SetAcceptance(m_etamin,m_etamax);
  return tracker;
}


void Tracker::Reset() {
  while (!m_tracks.empty()) { delete m_tracks.front(); m_tracks.pop_front(); }
  m_tracks.clear();
}

bool Tracker::Fill(const double E,const double eta,
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

void Tracker::Print(kf_code kfcode) {
  bool printit(true);
  for (std::list<Track *>::iterator trit=m_tracks.begin();
       trit!=m_tracks.end();trit++) {
    if (kfcode!=kf_none) {
      printit = false;
      if ((*trit)->flav.Kfcode()==kfcode) printit = true;
    }
    if (printit) {
      msg_Out()<<" Tracker : "<<(*trit)->flav<<" : "
	       <<(*trit)->mom[0]<<","<<(*trit)->eta<<","<<(*trit)->phi<<std::endl;
    }
  }
}

