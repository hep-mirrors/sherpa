#ifndef Primitive_Hadron_Calorimeter_H
#define Primitive_Hadron_Calorimeter_H

#include "Primitive_Detector.H"
#include "Particle_List.H"

namespace ATOOLS { class Particle_Qualifier_Base; }

namespace ANALYSIS {
  struct hadcalor {
    enum code {
      full     = 0,
      gaussian = 1
    };
  };

  struct hadgeom {
    enum code {
      full     = 0,
      flatmiss = 1
    };
  };

  class Primitive_Hadron_Calorimeter : public Primitive_Detector_Element {
  private:
    hadcalor::code m_calor;
    hadgeom::code  m_geom;
    double         m_thres, m_missprob, m_resolution;

    void   FillParticleInDetectorElement(const ATOOLS::Particle *);
    double EnergyDeposit(const ATOOLS::Particle *);
    void   SmearEnergy(const ATOOLS::Flavour & const,double &);
    void   GeometryEffects(const ATOOLS::Flavour & const,double &,
			   const double,const double);
    void   AddNoise();
  public:
    Primitive_Hadron_Calorimeter(const int,const int,
				 const double=-5.,const double=5.,
				 const double=0.,const double=2.*M_PI);
    ~Primitive_Hadron_Calorimeter();
    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    void Reset();
    void Print(std::ostream & =std::cout);

    void SetThreshold(const double thres)      { m_thres      = thres; }
    void SetMissProbabiltiy(const double miss) { m_missprob   = miss; }
    void SetEnergyResolution(const double res) { m_resolution = res; }

    void SetCalorimetryMode(const hadcalor::code hc) { m_calor = hc; }
    void SetGeometryMode(const hadgeom::code hg)     { m_geom  = hg; }

    Primitive_Detector_Element * Copy() const ;
  };
}

#endif
