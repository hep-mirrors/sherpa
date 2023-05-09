#ifndef PHOTONS_Photon_Splitter_Splitting_Functions_H
#define PHOTONS_Photon_Splitter_Splitting_Functions_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/PhotonSplitter/Kinematics_Base.H"
#include <cstddef>
#include <vector>

using namespace ATOOLS;

namespace PHOTONS {

  class Splitting_Function;
  class YFS_Particle;

  class YFS_Particle {
  private:
    size_t          m_id;
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_p;
    double          m_charge;

  public:
    YFS_Particle(const ATOOLS::Vec4D &p, const int &id, const ATOOLS::Flavour &flav);
    
    inline ATOOLS::Flavour GetFlavour() { return m_flav; }
    inline size_t Id() { return m_id; }
    inline ATOOLS::Vec4D Momentum() { return m_p; }
    inline void SetMomentum(ATOOLS::Vec4D newmom) { m_p = newmom; }
  };

  typedef std::vector<YFS_Particle *>                   YFS_Particle_Vector;
  typedef std::list<YFS_Particle *>                     YFS_Particle_List;

  class Splitting_Function {
  protected:
    size_t m_id;
    bool m_on;
    bool m_decay;
    YFS_Particle *p_splitter;
    double m_intspin;
    ATOOLS::Flavour m_flspec;
    ATOOLS::Flavour m_flavs[3];
    YFS_Particle_Vector m_specs;
    YFS_Particle * p_spec;
    Kinematics_FF      m_kinFF;
    Kinematics_FI      m_kinFI;
    Kinematics_IF      m_kinIF;
    Kinematics_II      m_kinII;
    double m_zmin, m_zmax;
    double m_mij2, m_mi2, m_mj2, m_mk2;
    double m_alpha, m_enhancefac;
    double m_startscale;
    double m_chargeCorrelator;

  public:
    Splitting_Function(YFS_Particle *splitter, int fla, 
        int flb, int flc, int intspin, const size_t &id, const double &enh, bool decay);

    ~Splitting_Function() {};

    void SetSpec(YFS_Particle *spec);
    void AddSpec(YFS_Particle *spec);
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
    inline YFS_Particle* GetSpectator() { return p_spec; }
    inline YFS_Particle_Vector GetSpecs() { return m_specs; }
    inline YFS_Particle* GetSplitter() { return p_splitter; }
    inline Kinematics_FF *KinFF() { return &m_kinFF; }
    inline Kinematics_FI *KinFI() { return &m_kinFI; }
    inline Kinematics_IF *KinIF() { return &m_kinIF; }
    inline Kinematics_II *KinII() { return &m_kinII; }

    double Lambda(const double &a, const double &b, const double &c) const;
    
    double JFF(const double Q2, const double y) const;

    double ChargeCorrelator();

    virtual double operator() (const double t, const double z, const double y, const double Q2) = 0;

    virtual double OverIntegrated(const double zmin,const double zmax) = 0;
    virtual double OverEstimated(const double z) = 0;
    virtual double Z() = 0;
  };

  class SF_FF : public Splitting_Function
  {
    public:
    SF_FF(YFS_Particle *splitter, int fla, 
        int flb, int flc, int intspin, const size_t &id, const double &enh,bool decay=1) :
        Splitting_Function(splitter,fla,flb,flc,intspin,id,enh,decay){};
    double operator()(const double t, const double z, const double y, const double Q2);
    double OverIntegrated(const double zmin,const double zmax);
    double OverEstimated(const double z);
    double Z();
  };

  class SF_FI : public Splitting_Function
  {
    public:
    SF_FI(YFS_Particle *splitter, int fla, 
        int flb, int flc, int intspin, const size_t &id, const double &enh,bool decay=1) :
        Splitting_Function(splitter,fla,flb,flc,intspin,id,enh,decay) {};
    double operator()(const double t, const double z, const double y, const double Q2);
    double OverIntegrated(const double zmin,const double zmax);
    double OverEstimated(const double z);
    double Z();
  };

  class SF_IF : public Splitting_Function
  {
    public:
    SF_IF(YFS_Particle *splitter, int fla, 
        int flb, int flc, int intspin, const size_t &id, const double &enh,bool decay=1) :
        Splitting_Function(splitter,fla,flb,flc,intspin,id,enh,decay) {};
    double operator()(const double t, const double z, const double y, const double Q2);
    double OverIntegrated(const double zmin,const double zmax);
    double OverEstimated(const double z);
    double Z();
  };

  class SF_II : public Splitting_Function
  {
    public:
    SF_II(YFS_Particle *splitter, int fla, 
        int flb, int flc, int intspin, const size_t &id, const double &enh,bool decay=1) :
        Splitting_Function(splitter,fla,flb,flc,intspin,id,enh,decay) {};
    double operator()(const double t, const double z, const double y, const double Q2);
    double OverIntegrated(const double zmin,const double zmax);
    double OverEstimated(const double z);
    double Z();
  };

  typedef std::vector<Splitting_Function *>             SF_Vector;
}

#endif
