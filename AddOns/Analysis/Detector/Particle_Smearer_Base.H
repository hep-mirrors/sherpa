#ifndef Analysis_Detector_Particle_Smearer_Base_H
#define Analysis_Detector_Particle_Smearer_Base_H

#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "ATOOLS/Phys/Particle_Qualifier.H"
#include "ATOOLS/Phys/Particle.H"
#include <map>

namespace ANALYSIS {
  class Primitive_Analysis;

  typedef std::pair<double,double> etainterval;

  class Particle_Smearer_Base : public Analysis_Object {
  protected:
    int    m_Emode_ECal, m_Dmode_ECal, m_Emode_HCal, m_Dmode_HCal;
    double m_missprob,m_threshold,m_evis_mean,m_evis_slope,m_efrac;
    double m_evis;

    ATOOLS::Particle * p_part;
    bool   m_track, m_mc;
    double m_eta_Track, m_phi_Track;
    double m_E_ECal, m_eta_ECal, m_phi_ECal;
    double m_E_HCal, m_eta_HCal, m_phi_HCal;
    double m_eta_MChamb, m_phi_MChamb;
    double m_E_deposed;

    std::string m_name;

    ATOOLS::Particle_Qualifier_Base * p_qualifier;

    std::map<etainterval, std::vector<double > *> m_Esmearingparams_ECal;
    std::map<etainterval, std::vector<double > *> m_Dsmearingparams_ECal;
    std::map<etainterval, std::vector<double > *> m_Esmearingparams_HCal;
    std::map<etainterval, std::vector<double > *> m_Dsmearingparams_HCal;
    std::map<etainterval, std::vector<double > *>::iterator m_spiter;

    inline bool IsInInterval(const double eta,const etainterval interval) { 
      return (eta>=interval.first && eta<=interval.second); 
    }
    void SetSmearingParams(etainterval,std::vector<double> *,
			   std::map<etainterval, std::vector<double > *> *);

    virtual void CalculateEvis(const double E);

    virtual void DetermineTracker()     =0;
    virtual void DetermineECal()        =0;
    virtual void DetermineHCal()        =0;
    virtual void DetermineMuonChambers()=0;


    virtual void Deflect(const bool,const int,const double,const double,const double);
    virtual void SmearEnergy(const bool,const int,const double,const double);

    std::vector<double > * GetSmearingParameters(std::map<etainterval, std::vector<double > *> *,
						 const double eta);
  public:
    Particle_Smearer_Base(Primitive_Analysis *,const std::string);
    virtual ~Particle_Smearer_Base();

    virtual void Evaluate(const ATOOLS::Blob_List &,double,double) {}
    virtual Analysis_Object *GetCopy() const { return NULL; }

    virtual bool   TreatParticle(ATOOLS::Particle *);
    virtual bool   GivesTrack(double & eta,double & phi);
    virtual double EnergyInECal(double & eta,double & phi);
    virtual double EnergyInHCal(double & eta,double & phi);
    virtual bool   GivesMuon(double & eta,double & phi);


    bool FillSegmentParameters(const bool mode,const bool cal,
			       std::vector<std::string> cur,const size_t);

    void SetEModes(const int m1,const int m2) { m_Emode_ECal = m1; m_Emode_HCal = m2; }
    void SetDModes(const int m1,const int m2) { m_Dmode_ECal = m1; m_Dmode_HCal = m2; }
    void SetMissProb(const double miss) { m_missprob = miss; }
    void SetThreshold(const double thr) { m_threshold = thr; }
    void SetEfrac(const double efrac)   { m_efrac = efrac; }
    void SetEvisfracParams(const double evis,const double eslope) {
      m_evis_mean = evis;m_evis_slope = eslope; 
    }
    

    const std::string Name() const { return m_name; }
  };
}

#endif
