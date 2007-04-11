#ifndef Primitive_Hadron_Calorimeter_H
#define Primitive_Hadron_Calorimeter_H

#include "Primitive_Detector_Element.H"
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
    double         m_threshold, m_missprob, m_resolution;

    bool   FillParticleInDetectorElement(const ATOOLS::Particle *);
    double EnergyDeposit(const ATOOLS::Particle *);
    void   SmearEnergy(const ATOOLS::Flavour &,double &);
    void   GeometryEffects(const ATOOLS::Flavour &,double &,
			   const double,const double);
    void   AddNoise();
  public:
    Primitive_Hadron_Calorimeter(const long int,const long int,
				 const double=-5.,const double=5.,
				 const double=0.,const double=2.*M_PI);
    ~Primitive_Hadron_Calorimeter();

    Analysis_Object * GetCopy() const;

    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    void Print(std::ostream & =std::cout);

    void SetThreshold(const double thres)      { m_threshold  = thres; }
    void SetMissProbability(const double miss) { m_missprob   = miss; }
    void SetResolution(const double res)       { m_resolution = res; }

    void SetCalorimetryMode(const hadcalor::code hc) { m_calor = hc; }
    void SetGeometryMode(const hadgeom::code hg)     { m_geom  = hg; }

    Primitive_Detector_Element * Copy() const ;
  };
}

#endif
