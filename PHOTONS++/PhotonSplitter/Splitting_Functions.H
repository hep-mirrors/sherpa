#ifndef PHOTONS_Photon_Splitter_Splitting_Functions_H
#define PHOTONS_Photon_Splitter_Splitting_Functions_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include <cstddef>
#include <vector>

using namespace ATOOLS;

namespace PHOTONS {

  class Splitting_Function;
  class Spectator;

  class Spectator {
  private:
    size_t          m_id;
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_p;
    double          m_charge;

  public:
    Spectator(const ATOOLS::Vec4D &p, const int &id, const ATOOLS::Flavour &flav);
    
    inline ATOOLS::Flavour GetFlavour() { return m_flav; }
    inline size_t Id() { return m_id; }
    inline ATOOLS::Vec4D Momentum() { return m_p; }
    inline void SetMomentum(ATOOLS::Vec4D newmom) { m_p = newmom; }
  };

  class Splitting_Function {
  private:
    size_t m_id;
    bool m_on;
    ATOOLS::Particle *p_splitter;
    double m_intspin;
    ATOOLS::Flavour m_flspec;
    ATOOLS::Flavour m_flavs[3];
    std::vector<Spectator*> m_specs;
    Spectator * p_spec;
    double m_zmin, m_zmax;
    double m_mij2, m_mi2, m_mj2, m_mk2;
    double m_alpha, m_enhancefac;
    double m_startscale;

  public:
    Splitting_Function(ATOOLS::Particle *splitter, int fla, 
        int flb, int flc, int intspin, const size_t &id, const double &enh);

    void SetSpec(Spectator *spec);
    void AddSpec(Spectator *spec);
    inline void SetLimits(const double zmin, const double zmax) { m_zmin = zmin; m_zmax = zmax; }
    inline double Cutoff() { return sqr(sqrt(m_mi2)+sqrt(m_mj2)); }
    inline double StartScale() { return m_startscale; }
    inline void SetStartScale(double tstart) { m_startscale = tstart; }
    inline bool On() { return m_on; }
    inline void SetOn(bool on) { m_on = on; }
    inline size_t Id() { return m_id; }
    inline ATOOLS::Flavour FlA() { return m_flavs[0]; }
    inline ATOOLS::Flavour FlB() { return m_flavs[1]; }
    inline ATOOLS::Flavour FlC() { return m_flavs[2]; }
    inline double Mass2A() { return m_mij2; }
    inline double Mass2B() { return m_mi2; }
    inline double Mass2C() { return m_mj2; }
    inline double Mass2Spec() { return m_mk2; }
    inline Spectator* GetSpectator() { return p_spec; }
    inline std::vector<Spectator*> GetSpecs() { return m_specs; }
    inline ATOOLS::Particle* GetSplitter() { return p_splitter; }

    double Lambda(const double &a, const double &b, const double &c) const;
    
    double JFF(const double Q2, const double y) const;
    
    double operator() (const double t, const double z, const double y, const double Q2);

    double OverIntegrated(const double zmin,const double zmax);
    double OverEstimated();
    double Z();
  };

  typedef std::vector<Splitting_Function *>   SF_Vector;
  typedef std::vector<Spectator *>            Spec_Vector;
}

#endif
