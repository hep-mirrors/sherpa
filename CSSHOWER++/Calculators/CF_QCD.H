#ifndef CSSHOWER_Main_CF_QCD_H
#define CSSHOWER_Main_CF_QCD_H

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"

namespace CSSHOWER { 
  const double s_Nc = 3.;
  const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
  const double s_CA = s_Nc;
  const double s_TR = 1./2.;

  class CF_QCD: public SF_Coupling {
  protected:
    class QCD_Coupling_Info {
    public:
      QCD_Coupling_Info(): p_cpl(NULL), m_rsf(1.0), p_cplmax(NULL) {}
      QCD_Coupling_Info(MODEL::Running_AlphaS * cpl, double rsf,
			std::vector<double> * cplmax=NULL):
	p_cpl(cpl), m_rsf(rsf), p_cplmax(cplmax) {}

      MODEL::Running_AlphaS * const Coupling() const { return p_cpl; }
      double RSF() const { return m_rsf; }
      std::vector<double> * const MaxCoupling() const { return p_cplmax; }
      void SetMaxCoupling(std::vector<double> * const cplmax) { p_cplmax = cplmax; }
      bool IsValid() const { return (p_cpl != NULL); }
    private:
      MODEL::Running_AlphaS * p_cpl;
      double m_rsf;
      std::vector<double> * p_cplmax;
    };

    /*!
     * Underlying couplings set by SetCoupling and
     * SetAlternativeUnderlyingCoupling, respectively
     *
     * If the alternative coupling is set, it takes precedence. This can
     * be used for reweighting purposes.
     */
    QCD_Coupling_Info m_maincplinfo, m_altcplinfo;
    const QCD_Coupling_Info & CurrentCouplingInfo() const
    { return (m_altcplinfo.IsValid() ? m_altcplinfo : m_maincplinfo); }

    //! Buffer of max alphas values to avoid re-calculations
    std::map<MODEL::Running_AlphaS *, std::vector<double> > m_altcplmax;
    double m_q, m_k0sq, m_kfac[7];
    int m_scvmode;
  public:
    CF_QCD(const SF_Key &key);
    double B0(const double &nf) const { return 11.0/6.0*s_CA-2.0/3.0*s_TR*nf; }
    bool SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    double CplMax(MODEL::Running_AlphaS * as, double rsf) const;
    double Coupling(const double &scale,const int pol);
    bool AllowSpec(const ATOOLS::Flavour &fl,const int mode);
    double CplFac(const double &scale) const;

    bool AllowsAlternativeCouplingUsage() const { return true; }
    void SetAlternativeUnderlyingCoupling(void *, double sf);
  };
}

#endif
