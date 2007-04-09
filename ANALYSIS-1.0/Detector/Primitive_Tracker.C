#include "Primitive_Tracker.H"
#include "Random.H"
#include "MathTools.H"

using namespace ANALYSIS;
using namespace ATOOLS;


Primitive_Tracker::
Primitive_Tracker(const int neta,const int nphi,
			     const double mineta,const double maxeta,
			     const double minphi,const double maxphi) :
  Primitive_Detector_Element(neta,nphi,mineta,maxeta,minphi,maxphi),
  m_trackprob(trackprob::full),m_threshold(0.),m_missprob(0.)
{
  p_qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject("charged","Charged"); 
  m_name      = "Tracker";
}


Primitive_Tracker::~Primitive_Tracker() {}

void Primitive_Tracker::Fill(const Particle_List * plist) {
  Reset();
  for (Particle_List::const_iterator pit=plist->begin();pit!=plist->end();pit++) {
    if (!(*p_qualifier)((*pit))) continue;
    FillParticleInDetectorElement((*pit));
  }
}

void Primitive_Tracker::
FillParticleInDetectorElement(const Particle * part) {
  int etapos,phipos;
  MatchCell(part->Momentum(),etapos,phipos);
  if (etapos<0||phipos<0) return;
  p_cells[etapos][phipos] = 
    TrackProbability(part->Flav(),Vec3D(part->Momentum()).Abs(),etapos,phipos);
}

double Primitive_Tracker::
TrackProbability(const Flavour & flav,const double p,const double eta,const double phi) {
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

void Primitive_Tracker::Reset()
{
  for (int i=0; i<m_neta; ++i) {
    for (int j=0; j<m_nphi; ++j) p_cells[i][j]=0.;
  }
  AddNoise();
}

void Primitive_Tracker::Print(std::ostream &) {}

void Primitive_Tracker::AddNoise() {}
