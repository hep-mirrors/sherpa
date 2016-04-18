#ifndef SHRiMPS_Soft_Jet_Criterion_H
#define SHRiMPS_Soft_Jet_Criterion_H

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PDF/Main/Jet_Criterion.H"
#include <map>

namespace SHRIMPS {
  class Soft_Jet_Criterion : public PDF::Jet_Criterion {
    ATOOLS::Cluster_Amplitude * p_refampl;
    ATOOLS::ClusterLeg_Vector   m_reflegs; 
    std::map<ATOOLS::Cluster_Leg *,double> m_kt2veto;

    size_t m_orig;
    double m_pt2;
    
    size_t FindCombination(ATOOLS::Cluster_Amplitude *ampl);
    double PTij2(ATOOLS::Vec4D pi,ATOOLS::Vec4D pj);
  public:
    Soft_Jet_Criterion(): p_refampl(NULL) {}
    ~Soft_Jet_Criterion() {}

    inline void SetClusterAmplitude(ATOOLS::Cluster_Amplitude * ampl) {
      p_refampl = ampl;
      m_reflegs = p_refampl->Legs();
    } 
    inline void Reset() {
      m_kt2veto.clear();
      m_reflegs.clear();
      p_refampl = NULL;
    }
    inline void SetKT2Veto(ATOOLS::Cluster_Leg * leg,double kt2veto) {
      m_kt2veto[leg] = kt2veto;
    }

    bool Jets(ATOOLS::Cluster_Amplitude *ampl,int mode=0);
    void Output();
  };


  class JF : public PHASIC::Jet_Finder {
  private:
    double m_ycut;
  public:
    JF();
    ~JF() {}

    double Ycut()             { return m_ycut; }
    void   SetJetCriterion(PDF::Jet_Criterion * jc) { p_jc = jc; }
  };
}

#endif
