#ifndef MCATNLO_Showers_SF_Coupling_H
#define MCATNLO_Showers_SF_Coupling_H

#include "MCATNLO/Showers/SF_Key.H"
#include "MCATNLO/Tools/Parton.H"

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {
  class Cluster_Amplitude;
}

namespace MCATNLO {

  class SF_Lorentz;

  class SF_Coupling {
  protected:

    SF_Lorentz *p_lf;

    cstp::code m_type;
    double     m_cplfac;

    //! Last unpolarized coupling value
    double     m_last;

    int        m_kfmode;

    /**
      * m_cplmax[i] = max_{scale} Coupling(scale, i)
      */
    std::vector<double> m_cplmax;

    static double s_qfac;

  public:

    SF_Coupling(const SF_Key &key);

    virtual ~SF_Coupling();

    virtual bool SetCoupling(MODEL::Model_Base *md,
			     const double &k0sqi,const double &k0sqf,
			     const double &isfac,const double &fsfac) = 0;

    //! Subclasses can return true if they support alternative couplings to be used (e.g. for reweighting)
    virtual bool AllowsAlternativeCouplingUsage() const { return false; }
    /*!
     * Use an alternative coupling
     *
     * @param cpl See subclass implementation(s) for the accepted pointer type.
     * @param sf An optional scale factor, see subclass implementations for the
     * concrete meaning of this argument.
     * Passing NULL ends the use of the alternative coupling.
     */
    virtual void SetAlternativeUnderlyingCoupling(void * cpl, double sf=1.0) {}

    virtual double Coupling(const double &scale,const int pol,
			    ATOOLS::Cluster_Amplitude *const sub) = 0;
    virtual bool AllowSpec(const ATOOLS::Flavour &fl) = 0;

    //! Last unpolarized coupling value
    double Last() const { return m_last; }

    //! Set last unpolarized coupling value
    void SetLast(double last) { m_last = last; }

    virtual double CplFac(const double &scale) const;

    virtual void ColorPoint(Parton *const p) const;

    virtual double ColorWeight(const Color_Info &ci) const;

    inline double MaxCoupling(const int mode) const { return s_qfac*m_cplmax[mode]; }

    inline void SetLF(SF_Lorentz *const lf) { p_lf=lf; }

    inline int KFMode() { return m_kfmode; }

    inline void SetKFMode(const int kfmode) { m_kfmode=kfmode; }

    inline static void SetQFac(const double &qfac) { s_qfac=qfac; }

  };

  typedef ATOOLS::Getter_Function<SF_Coupling,SF_Key,
				  std::less<std::string> > SFC_Getter;

  // NOTE: SFC_Filler_Key/Getter allow for specifying SF_Coupling Getters
  // (SFC_Getter) during run-time; this is triggered via calling
  // SVC_Filler_Getter's GetObject in Sudakov::InitSplittingFunctions.
  // This usage is highly unintuitive (and should probably be refactored),
  // therefore we document it here.

  struct SFC_Filler_Key {
    const MODEL::Model_Base *p_md;
    std::vector<SFC_Getter*> *p_gets;
    inline SFC_Filler_Key(const MODEL::Model_Base *const md,
			  std::vector<SFC_Getter*> *const gets):
      p_md(md), p_gets(gets) {}
  };

  typedef ATOOLS::Getter_Function
  <void,SFC_Filler_Key,std::less<std::string> > SFC_Filler_Getter;

}

#define DECLARE_CPL_GETTER(NAME)					\
									\
  class NAME: public ATOOLS::Getter_Function				\
  <SF_Coupling,SF_Key,std::less<std::string> > {			\
  protected:								\
    void PrintInfo(std::ostream &str,const size_t width) const;		\
    Object_Type *							\
      operator()(const Parameter_Type &parameters) const;		\
  public:								\
    NAME(const std::string &name):					\
      ATOOLS::Getter_Function<SF_Coupling,SF_Key,			\
      std::less<std::string> >(name) {}					\
  }

#endif
