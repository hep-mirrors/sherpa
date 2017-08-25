#ifndef PDF__Main__Cluster_Definitions_Base_H
#define PDF__Main__Cluster_Definitions_Base_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Decay_Info.H"
#include "ATOOLS/Math/Poincare.H"

namespace PDF {

  struct Cluster_Config {
    ATOOLS::Cluster_Amplitude *p_ampl;
    int m_i, m_j, m_k, m_kin, m_mode;
    ATOOLS::Flavour m_mo;
    const ATOOLS::Mass_Selector *p_ms;
    inline Cluster_Config
    (ATOOLS::Cluster_Amplitude *ampl,
     const int i,const int j=-1,const int k=-1,
     const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
     const ATOOLS::Mass_Selector *ms=NULL,
     ATOOLS::Decay_Info *const dec=NULL,
     const int kin=-1,const int mode=0):
      p_ampl(ampl), m_i(i), m_j(j), m_k(k),
      m_kin(kin), m_mode(mode),
      m_mo(mo), p_ms(ms) {}
    bool PureQCD() const;
    bool operator<(const Cluster_Config &cc) const;
  };// end of struct Cluster_Config

  std::ostream &operator<<(std::ostream &str,const Cluster_Config &cc);

  class Cluster_Definitions_Base;

  struct Cluster_Param {
    Cluster_Definitions_Base *p_ca;
    double m_op, m_kt2, m_mu2;
    int m_cpl, m_kin, m_mode, m_stat;
    ATOOLS::Vec4D m_pijt, m_pkt;
    ATOOLS::Poincare_Sequence m_lam;
    Cluster_Param(Cluster_Definitions_Base *const ca=NULL,
		  const double &op=0.0,
		  const double &kt2=0.0,const double &mu2=0.0,
		  const int cpl=-1,const int kin=0,const int mode=0,
		  const ATOOLS::Vec4D &pijt=ATOOLS::Vec4D(),
		  const ATOOLS::Vec4D &pkt=ATOOLS::Vec4D(),
		  const ATOOLS::Poincare_Sequence &lam=
		  ATOOLS::Poincare_Sequence()):
      p_ca(ca), m_op(op), m_kt2(kt2), m_mu2(mu2),
      m_cpl(cpl), m_kin(kin), m_mode(mode), m_stat(0),
      m_pijt(pijt), m_pkt(pkt), m_lam(lam) {}
  };// end of struct Cluster_Param

  std::ostream &operator<<(std::ostream &str,const Cluster_Param &cp);

  typedef std::pair<Cluster_Config,Cluster_Param> Cluster_Info;

  typedef std::vector<Cluster_Info> ClusterInfo_Vector;

  class Cluster_Definitions_Base {
  public:

    Cluster_Definitions_Base();
    virtual ~Cluster_Definitions_Base();

    virtual Cluster_Param Cluster(const Cluster_Config &ca) = 0;

    int ReCluster(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Cluster_Definitions_Base

}// end of namespace PDF

#endif
