#include "Primitive_Tracker.H"
#include "Random.H"
#include "MathTools.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Primitive_Tracker_Getter,"Tracker",Analysis_Object,Argument_Matrix);

Analysis_Object *				     
Primitive_Tracker_Getter::operator()(const Argument_Matrix & params) const  
{ 
  double   mineta(-5.), maxeta(5.), minphi(0.), maxphi(2.*M_PI);
  long int neta(1000000), nphi(628319);
  trackprob::code tp=trackprob::full;
  double threshold(0.),missprob(0.);
  for (size_t i=0;i<params.size();++i) {
    const std::vector<std::string> &cur=params[i];
    if (cur[0]=="Acceptance" && cur.size()==5) {
      mineta=ToType<double>(cur[1]);maxeta=ToType<double>(cur[2]);
      minphi=ToType<double>(cur[3]);maxphi=ToType<double>(cur[4]);
    }
    else if (cur[0]=="Granularity" && cur.size()==3) {
      neta=ToType<long int>(cur[1]);nphi=ToType<long int>(cur[2]);
    }
    else if (cur[0]=="Tracking_Mode" && cur.size()==2) {
      tp=trackprob::code(ToType<int>(cur[1]));
    }
    else if (cur[0]=="Miss_Probability" && cur.size()==2 &&
	     (tp==trackprob::flatmiss || tp==trackprob::threshold_flatmiss)) {
      missprob=ToType<double>(cur[1]);
    }
    else if (cur[0]=="Tracking_Threshold" && cur.size()==2 &&
	     (tp==trackprob::threshold || tp==trackprob::threshold_flatmiss)) {
      threshold=ToType<double>(cur[1]);
    }
  }
  Primitive_Tracker * tracker(new Primitive_Tracker(neta,nphi,mineta,maxeta,minphi,maxphi));
  tracker->SetTrackingMode(tp);
  tracker->SetThreshold(threshold);
  tracker->SetMissProbability(missprob);
  return tracker; 
}

void Primitive_Tracker_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Tracker"; 
}

Primitive_Tracker::
Primitive_Tracker(const long int neta,const long int nphi,
		  const double mineta,const double maxeta,
		  const double minphi,const double maxphi) :
  Primitive_Detector_Element(pde::tracks,neta,nphi,mineta,maxeta,minphi,maxphi),
  m_trackprob(trackprob::full),m_threshold(0.),m_missprob(0.)
{
  p_qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject("charged","Charged"); 
  m_name      = "Tracker";
}


Primitive_Tracker::~Primitive_Tracker() {}

Analysis_Object * Primitive_Tracker::GetCopy() const {
  Primitive_Tracker * tracker =
    new Primitive_Tracker(m_neta,m_nphi,m_etamin,m_etamax,m_phimin,m_phimax);
  tracker->SetMissProbability(m_missprob);
  tracker->SetThreshold(m_threshold);
  tracker->SetTrackingMode(m_trackprob);
  return tracker;
}

void Primitive_Tracker::Fill(const Particle_List * plist) {
  Reset();
  for (Particle_List::const_iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!(*p_qualifier)((*pit))) continue;
    FillParticleInDetectorElement((*pit));
  }
}

bool Primitive_Tracker::FillParticleInDetectorElement(const Particle * part) {
  long int etapos,phipos;
  MatchCell(part->Momentum(),etapos,phipos);
  if (etapos<0||phipos<0) return false;
  std::pair<long int, long int> pos;
  pos.first=etapos; pos.second=phipos;
  double trackprob =
    TrackProbability(part->Flav(),Vec3D(part->Momentum()).Abs(),etapos,phipos);
  if (m_tracks.find(pos)!=m_tracks.end()) 
    m_tracks.find(pos)->second += trackprob;
  else m_tracks[pos] = trackprob;
  return false;
}

double Primitive_Tracker::TrackProbability(const Flavour & flav,const double p,
					   const double eta,const double phi) {
  switch (m_trackprob) {
    case trackprob::threshold_flatmiss:
      if (p>m_threshold && ran.Get()<m_missprob) return 1.;
      break;
    case trackprob::threshold:
      if (p>m_threshold) return 1.;
      break;
    case trackprob::flatmiss:
      if (ran.Get()<m_missprob) return 1.;
      break;
    case trackprob::full:
    default:
      return 1.;
  }
  return 0.;
}

void Primitive_Tracker::Extract(Particle_List *) {}

void Primitive_Tracker::Print(std::ostream &) {}

void Primitive_Tracker::AddNoise() {}
