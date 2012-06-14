#ifndef POWHEG_Main_CS_Gamma_H
#define POWHEG_Main_CS_Gamma_H

#include "POWHEG/Showers/Shower.H"
#include "POWHEG/Main/CS_Cluster_Definitions.H"
#include "PDF/Main/NLOMC_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace PHASIC {
  class Process_Base;
  class Single_Process;
}

namespace POWHEG {

  struct Weight_Key {
    size_t m_ij, m_k;
    ATOOLS::RB_Key m_rbkey;
    Weight_Key(const size_t &ij,const size_t &k,const ATOOLS::Flavour &flij,
	       const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj);
    inline bool operator<(const Weight_Key &wk) const
    { return m_ij<wk.m_ij?true:(m_ij>wk.m_ij?false:m_k<wk.m_k); }
  };// end of struct Weight_Key

  std::ostream &operator<<(std::ostream &str,const Weight_Key &k);

  struct Weight_Value {
    PHASIC::Process_Base *p_proc;
    Splitting_Function_Base *p_sf;
    double m_me, m_b, m_qij2, m_muf2, m_mur2;
    bool m_trig;
    inline Weight_Value(PHASIC::Process_Base *const proc=NULL,
			const double &muf2=0.0,const double &mur2=0.0):
      p_proc(proc), p_sf(NULL), m_me(0.0), m_b(0.0), m_qij2(0.0),
      m_muf2(muf2), m_mur2(mur2), m_trig(true) {}
  };// end of struct Weight_Value

  std::ostream &operator<<(std::ostream &str,const Weight_Value &w);

  typedef std::map<Weight_Key,Weight_Value> Weight_Map;

  class CS_POWHEG;

  class CS_Gamma {
  private:

    CS_POWHEG *p_css;
    Shower    *p_shower;

    CS_Cluster_Definitions *p_cluster;
    ATOOLS::Mass_Selector  *p_ms;

    PHASIC::Single_Process *p_rproc;

    int    m_on;
    double m_zhth, m_ktres;

    bool CheckCore(const ATOOLS::Cluster_Amplitude *ampl) const;

    bool CheckColors(const ATOOLS::Cluster_Leg *li,
		     const ATOOLS::Cluster_Leg *lj,
		     const ATOOLS::Cluster_Leg *lk,
		     const ATOOLS::Flavour &mo) const;

    int SingleWeight(ATOOLS::Cluster_Amplitude *const rampl,
		     ATOOLS::Cluster_Leg *const li,
		     ATOOLS::Cluster_Leg *const lj,
		     ATOOLS::Cluster_Leg *const lk,const CS_Parameters &cs,
		     const std::map<size_t,size_t> &idmap,Weight_Map &ws,
		     const int mode,const bool userb);
    int CalculateWeights(ATOOLS::Cluster_Amplitude *const rampl,
			 const std::map<size_t,size_t> &idmap,
			 Weight_Map &ws,const int mode,const bool userb);

    Weight_Map CalculateWeight(ATOOLS::Cluster_Amplitude *const ampl,
			       const int mode=0,const bool userb=true);

    Weight_Value Differential(ATOOLS::Cluster_Amplitude *const ampl,
			      const int mode=2) const;

    double TrialWeight(ATOOLS::Cluster_Amplitude *const ampl);

  public:

    CS_Gamma(CS_POWHEG *const css,Shower *const shower,
	     CS_Cluster_Definitions *const clus);

    bool Reject();
    void AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl);

    inline void SetZHParams(const double &zhth,const double &ktres)
    { m_zhth=zhth; m_ktres=ktres; }

    inline void SetOn(const int on) { m_on=on; }

    inline int On() const { return m_on; }

  };// end of class CS_Gamma

}// end of namespace POWHEG

#endif
