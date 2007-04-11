#ifndef Primitive_ElMag_Calorimeter_H
#define Primitive_ElMag_Calorimeter_H

#include "Primitive_Detector_Element.H"
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
    std::vector<std::vector<int> > m_id;

    emcalor::code m_calor;
    emgeom::code  m_geom;
    double        m_threshold, m_missprob, m_resolution;

    ATOOLS::Particle_Qualifier_Base *p_electron, *p_photon, *p_muon, *p_charged_hadron;

    bool   FillParticleInDetectorElement(const ATOOLS::Particle *);
    double EnergyDeposit(const ATOOLS::Particle *,bool &);
    void   SmearEnergy(const ATOOLS::Flavour &,double &);
    void   GeometryEffects(const ATOOLS::Flavour &,double &,
			   const double,const double);
    void   AddNoise();
  public:
    Primitive_ElMag_Calorimeter(const long int,const long int,
				const double=-5.,const double=5.,
				const double=0.,const double=2.*M_PI);
    ~Primitive_ElMag_Calorimeter();

    Analysis_Object * GetCopy() const;

    void Reset();
    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    void Print(std::ostream & =std::cout);
    int  GetFlav(const long int, const long int) const;

    void SetThreshold(const double thres)      { m_threshold  = thres; }
    void SetMissProbability(const double miss) { m_missprob   = miss; }
    void SetResolution(const double res)       { m_resolution = res; }

    void SetCalorimetryMode(const emcalor::code elc) { m_calor = elc; }
    void SetGeometryMode(const emgeom::code elg)     { m_geom  = elg; }

    Primitive_Detector_Element * Copy() const ;
  };
}

#endif
