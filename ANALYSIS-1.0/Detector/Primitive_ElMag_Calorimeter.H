#ifndef Primitive_ElMag_Calorimeter_H
#define Primitive_ElMag_Calorimeter_H

#include "Primitive_Detector.H"
#include "Particle_List.H"

namespace ATOOLS { class Particle_Qualifier_Base; }

namespace ANALYSIS {
  struct emcalor {
    enum code {
      full     = 0,
      gaussian = 1
    };
  };

  struct emgeom {
    enum code {
      full     = 0,
      flatmiss = 1
    };
  };

  class Primitive_ElMag_Calorimeter : public Primitive_Detector_Element {
  private:
    emcalor::code m_calor;
    emgeom::code  m_geom;
    double        m_thres, m_missprob, m_resolution;

    ATOOLS::Particle_Qualifier * p_electron, * p_photon, * p_muon, * p_charged_hadron;

    void   FillParticleInDetectorElement(const ATOOLS::Particle *);
    double EnergyDeposit(const ATOOLS::Particle *);
    void   SmearEnergy(const ATOOLS::Flavour & const,double &);
    void   GeometryEffects(const ATOOLS::Flavour & const,double &,
			   const double,const double);
    void   AddNoise();
  public:
    Primitive_ElMag_Calorimeter(const int,const int,
				const double=-5.,const double=5.,
				const double=0.,const double=2.*M_PI);
    ~Primitive_ElMag_Calorimeter();
    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    void Reset();
    void Print(std::ostream & =std::cout);

    void SetThreshold(const double thres)      { m_thres      = thres; }
    void SetMissProbabiltiy(const double miss) { m_missprob   = miss; }
    void SetEnergyResolution(const double res) { m_resolution = res; }

    void SetCalorimetryMode(const emcalor::code elc) { m_calor = elc; }
    void SetGeometryMode(const emgeom::code elg)     { m_geom  = elg; }

    Primitive_Detector_Element * Copy() const ;
  };
}

#endif
