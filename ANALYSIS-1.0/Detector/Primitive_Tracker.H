#ifndef Primitive_Tracker_H
#define Primitive_Tracker_H

#include "Primitive_Detector_Element.H"
#include "Particle_List.H"

namespace ATOOLS { class Particle_Qualifier_Base; }

namespace ANALYSIS {
  struct trackprob {
    enum code {
      full               = 0,
      flatmiss           = 1,
      threshold          = 2,
      threshold_flatmiss = 3
    };
  };

  class Primitive_Tracker : public Primitive_Detector_Element {
  private:
    trackprob::code m_trackprob;
    double          m_threshold, m_missprob;

    void   FillParticleInDetectorElement(const ATOOLS::Particle *);
    double TrackProbability(const ATOOLS::Flavour &,const double,
			    const double,const double);
    void   AddNoise();
  public:
    Primitive_Tracker(const int,const int,
		      const double=-5.,const double=5.,
		      const double=0.,const double=2.*M_PI);
    ~Primitive_Tracker();
    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    void Reset();
    void Print(std::ostream & =std::cout);

    void SetMissProbabilty(const double miss)     { m_missprob   = miss; }
    void SetTrackingMode(const trackprob::code tp) { m_trackprob  = tp; }

    Primitive_Detector_Element * Copy() const ;
  };
}

#endif
