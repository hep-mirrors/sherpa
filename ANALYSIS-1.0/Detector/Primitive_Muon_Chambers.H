#ifndef Primitive_Muon_Chambers_H
#define Primitive_Muon_Chambers_H

#include "Primitive_Detector_Element.H"
#include "Particle_List.H"

namespace ATOOLS { class Particle_Qualifier_Base; }

namespace ANALYSIS {
  struct muonprob {
    enum code {
      full     = 0,
      flatmiss = 1
    };
  };

  class Primitive_Muon_Chambers : public Primitive_Detector_Element {
  private:
    muonprob::code m_muonprob;
    double         m_threshold,m_missprob,m_resolution;

    void   FillParticleInDetectorElement(const ATOOLS::Particle *);
    double EnergyDeposit(const ATOOLS::Particle *);
    void   SmearEnergy(const ATOOLS::Flavour &,double &);
    void   GeometryEffects(const ATOOLS::Flavour &,double &,
			   const double,const double);
    void   AddNoise();
  public:
    Primitive_Muon_Chambers(const long int,const long int,
			    const double=-5.,const double=5.,
			    const double=0.,const double=2.*M_PI);
    ~Primitive_Muon_Chambers();
    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    void Print(std::ostream & =std::cout);

    void SetThreshold(const double thres)      { m_threshold  = thres; }
    void SetMissProbability(const double miss) { m_missprob   = miss; }
    void SetResolution(const double res)       { m_resolution = res; }

    void SetMuonMode(const muonprob::code mp)  { m_muonprob = mp; }

    Primitive_Detector_Element * Copy() const ;
  };
}

#endif
