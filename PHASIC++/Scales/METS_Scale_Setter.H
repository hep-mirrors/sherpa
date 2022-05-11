#ifndef PHASIC_Scales_METS_Scale_Setter_H
#define PHASIC_Scales_METS_Scale_Setter_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PHASIC++/Scales/Tag_Setter.H"

#include <limits>

namespace PHASIC {

  struct CS_Params {
    size_t m_idi, m_idj, m_idk, m_i, m_j, m_k, m_oqcd;
    ATOOLS::Flavour m_fl;
    double m_kt2, m_op2, m_mu2, m_z, m_y;
    ATOOLS::Decay_Info *p_dec;
    ATOOLS::Vec4D m_pijt, m_pkt;
    ATOOLS::Poincare_Sequence m_lam;
    CS_Params(const size_t &idi,const size_t &idj,
	      const size_t &idk,const size_t &i,const size_t &j,
	      const size_t &k,const ATOOLS::Flavour &fl):
      m_idi(idi),m_idj(idj), m_idk(idk),
      m_i(i),m_j(j), m_k(k), m_oqcd(0), m_fl(fl),
      m_kt2(-1.0), m_op2(-std::numeric_limits<double>::max()),
      m_mu2(-1.0), m_z(0.0), m_y(0.0), p_dec(NULL) {}
    bool operator<(const CS_Params &ck) const
    { 
      if (m_idi<ck.m_idi) return true;
      if (m_idi>ck.m_idi) return false;
      if (m_idj<ck.m_idj) return true;
      if (m_idj>ck.m_idj) return false;
      if (m_idk<ck.m_idk) return true;
      if (m_idk>ck.m_idk) return false;
      return m_fl<ck.m_fl;
    }
    void SetParams(const double &kt2,const double &z,const double &y,
		   const ATOOLS::Vec4D &pijt,const ATOOLS::Vec4D &pkt,
		   const ATOOLS::Poincare_Sequence &lam=
		   ATOOLS::Poincare_Sequence())
    { m_mu2=m_kt2=kt2, m_z=z; m_y=y; m_pijt=pijt; m_pkt=pkt; m_lam=lam; }
  };// end of struct CS_Params

  class METS_Scale_Setter: public Scale_Setter_Base {
  private:

    Core_Scale_Setter *p_core, *p_uoscale;

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

    ATOOLS::Flavour_Vector m_f;

    SP(Color_Integrator) p_ci;

    size_t m_cnt, m_rej, m_mode, m_cmode, m_cmodebvi, m_cmoders;
    double m_lfrac, m_aqed, m_wthres, m_rsf, m_csf;
    int    m_rproc, m_sproc, m_vproc, m_nproc, m_nfgsplit, m_allowuo;

    ATOOLS::DecayInfo_Vector m_decids;

    static double s_eps, s_kt2max;

    PDF::CParam UnorderedScale(ATOOLS::Cluster_Amplitude *const ampl) const;

    bool CheckColors(const ATOOLS::Cluster_Leg *li,
		     const ATOOLS::Cluster_Leg *lj,
		     const ATOOLS::Cluster_Leg *lk,
		     const ATOOLS::Flavour &mo) const;
    ATOOLS::ColorID CombineColors(const ATOOLS::Cluster_Leg *li,
				  const ATOOLS::Cluster_Leg *lj,
				  const ATOOLS::Cluster_Leg *lk,
				  const ATOOLS::Flavour &mo) const;

    double Lam(const double &s,
	       const double &sb,const double &sc) const;

    void KT2(const ATOOLS::Cluster_Leg *li,
	     const ATOOLS::Cluster_Leg *lj,
	     const ATOOLS::Cluster_Leg *lk,CS_Params &cs) const;

    bool Combine(ATOOLS::Cluster_Amplitude &ampl,
		 int i,int j,int k,const CS_Params &cs) const;

    double SetScales(const double &muf2,ATOOLS::Cluster_Amplitude *ampl,
		     const size_t &mode);

    double CalculateStrict(const ATOOLS::Vec4D_Vector &momenta,
			   const size_t &mode);

  public:

    METS_Scale_Setter(const Scale_Setter_Arguments &args,
		      const int mode=1);

    ~METS_Scale_Setter();

    PDF::CParam CoreScale(ATOOLS::Cluster_Amplitude *const ampl) const override;

    void PreCalc(const ATOOLS::Vec4D_Vector &p,const size_t &mode) override;

    double Calculate(const ATOOLS::Vec4D_Vector &p,const size_t &mode) override;

    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class METS_Scale_Setter

  class Loose_METS_Scale_Setter: public METS_Scale_Setter {};
  class Strict_METS_Scale_Setter: public METS_Scale_Setter {};

}// end of namespace PHASIC

#endif
