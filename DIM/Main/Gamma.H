#ifndef DIM__Main__Gamma_H
#define DIM__Main__Gamma_H

#include "DIM/Shower/Shower.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace PHASIC {
  class Process_Base;
  class Single_Process;
}

namespace DIM {

  struct Weight_Key {
    size_t m_ij, m_k;
    inline Weight_Key(const size_t &ij,const size_t &k):
      m_ij(ij), m_k(k) {}
    inline bool operator<(const Weight_Key &wk) const
    { return m_ij<wk.m_ij?true:(m_ij>wk.m_ij?false:m_k<wk.m_k); }
  };// end of struct Weight_Key

  std::ostream &operator<<(std::ostream &str,const Weight_Key &k);

  struct Weight_Value {
    PHASIC::Process_Base *p_proc;
    const Kernel *p_sf;
    double m_me, m_b, m_muf2[2], m_mur2, m_muq2;
    inline Weight_Value(PHASIC::Process_Base *const proc=NULL):
      p_proc(proc), p_sf(NULL), m_me(0.0), m_b(0.0),
      m_muf2{0.0,0.0}, m_mur2(0.0), m_muq2(0.0) {}
  };// end of struct Weight_Value

  std::ostream &operator<<(std::ostream &str,const Weight_Value &w);

  typedef std::map<Weight_Key,Weight_Value> Weight_Map;

  class MCatNLO;

  class Gamma {
  private:

    MCatNLO *p_dire;
    Shower  *p_shower;

    ATOOLS::Mass_Selector *p_ms;

    double m_weight;

    Weight_Value Differential(ATOOLS::Cluster_Amplitude *const ampl,
			      const ATOOLS::nlo_type::code type=ATOOLS::nlo_type::lo,
			      const std::string add="") const;

    Weight_Map CalculateWeight(ATOOLS::Cluster_Amplitude *const ampl);

    MC_Weight TrialWeight(ATOOLS::Cluster_Amplitude *const ampl);

  public:

    Gamma(MCatNLO *const dire,Shower *const shower);

    bool Reject();

    inline double Weight() const { return m_weight; }

  };// end of class Gamma

}// end of namespace DIM

#endif
