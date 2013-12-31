#ifndef MCATNLO_Main_CS_Gamma_H
#define MCATNLO_Main_CS_Gamma_H

#include "MCATNLO/Showers/Shower.H"
#include "MCATNLO/Main/CS_Cluster_Definitions.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace PHASIC {
  class Process_Base;
  class Single_Process;
}

namespace MCATNLO {

  struct Weight_Key {
    size_t m_ij, m_k;
    Weight_Key(const size_t &ij,const size_t &k,const ATOOLS::Flavour &flij,
	       const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj);
    inline bool operator<(const Weight_Key &wk) const
    { return m_ij<wk.m_ij?true:(m_ij>wk.m_ij?false:m_k<wk.m_k); }
  };// end of struct Weight_Key

  std::ostream &operator<<(std::ostream &str,const Weight_Key &k);

  struct Weight_Value {
    PHASIC::Process_Base *p_proc;
    Splitting_Function_Base *p_sf;
    double m_me, m_b, m_muf2, m_mur2;
    inline Weight_Value(PHASIC::Process_Base *const proc=NULL,
			const double &muf2=0.0,const double &mur2=0.0):
      p_proc(proc), p_sf(NULL), m_me(0.0), m_b(0.0),
      m_muf2(muf2), m_mur2(mur2) {}
  };// end of struct Weight_Value

  std::ostream &operator<<(std::ostream &str,const Weight_Value &w);

  typedef std::map<Weight_Key,Weight_Value> Weight_Map;

  class CS_MCatNLO;

  struct Trial_Weight {
  private:
    double m_f, m_g, m_h;
  public:
    inline Trial_Weight(const double &f,
			const double &g,const double &h):
      m_f(f), m_g(g), m_h(h) {}
    inline double MC() const { return m_f/m_g; }
    inline double Accept() const { return m_g/m_h; }
    inline double Reject() const
    { return m_g/m_h*(m_h-m_f)/(m_g-m_f); }
  };// end of struct Trial_Weight

  class CS_Gamma {
  private:

    CS_MCatNLO *p_css;
    Shower    *p_shower;

    CS_Cluster_Definitions *p_cluster;
    ATOOLS::Mass_Selector  *p_ms;

    int    m_on;
    double m_weight, m_oef;

    bool CheckColors(const ATOOLS::Cluster_Leg *li,
		     const ATOOLS::Cluster_Leg *lj,
		     const ATOOLS::Cluster_Leg *lk,
		     const ATOOLS::Flavour &mo) const;
    ATOOLS::ColorID CombineColors(const ATOOLS::Cluster_Leg *li,
				  const ATOOLS::Cluster_Leg *lj,
				  const ATOOLS::Cluster_Leg *lk,
				  const ATOOLS::Flavour &mo) const;

    int SingleWeight(ATOOLS::Cluster_Amplitude *const rampl,
		     ATOOLS::Cluster_Leg *const li,
		     ATOOLS::Cluster_Leg *const lj,
		     ATOOLS::Cluster_Leg *const lk,const CS_Parameters &cs,
		     const std::map<size_t,size_t> &idmap,Weight_Map &ws,
		     const int mode);
    int CalculateWeights(ATOOLS::Cluster_Amplitude *const rampl,
			 const std::map<size_t,size_t> &idmap,
			 Weight_Map &ws,const int mode);

    Weight_Map CalculateWeight(ATOOLS::Cluster_Amplitude *const ampl,
			       const int mode);

    Weight_Value Differential(ATOOLS::Cluster_Amplitude *const ampl,
			      const PHASIC::nlo_type::code type=PHASIC::nlo_type::lo,
			      const std::string add="") const;

    Trial_Weight TrialWeight(ATOOLS::Cluster_Amplitude *const ampl);

  public:

    CS_Gamma(CS_MCatNLO *const css,Shower *const shower,
	     CS_Cluster_Definitions *const clus);

    bool Reject();
    void AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl);

    inline void SetOn(const int on) { m_on=on; }

    inline int On() const { return m_on; }

    inline double Weight() const { return m_weight; }

    inline void SetOEF(const double &oef) { m_oef=oef; }

  };// end of class CS_Gamma

}// end of namespace MCATNLO

#endif
