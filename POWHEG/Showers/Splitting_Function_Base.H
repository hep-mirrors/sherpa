#ifndef POWHEG_Showers_Splitting_Function_Base_H
#define POWHEG_Showers_Splitting_Function_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/RB_Data.H"
#include "ATOOLS/Math/Vector.H"

namespace ATOOLS { 
  class Function_Base;
  class Cluster_Amplitude;
}

namespace MODEL {
  class Model_Base; 
  class Single_Vertex;
}

namespace PDF { class PDF_Base; }

namespace POWHEG {

  class Parton;

  struct cstp {
    enum code {
      FF   = 11,
      FI   = 12,
      IF   = 21,
      II   = 22,
      none = 0
    };
  };

  class SF_Coupling;

  struct SF_Key {
  public:
    ATOOLS::Mass_Selector *&p_ms;
    const MODEL::Single_Vertex *p_v;
    int m_mode, m_kfmode, m_ewmode;
    cstp::code m_type;
    SF_Coupling *p_cf;
  public:
    inline SF_Key(ATOOLS::Mass_Selector *&ms,
		  const MODEL::Single_Vertex *v,
		  const int mode,const cstp::code type,
		  const int kfmode,const int ewmode):
      p_ms(ms), p_v(v), m_mode(mode), m_kfmode(kfmode),
      m_ewmode(ewmode), m_type(type), p_cf(NULL) {}
    std::string ID(const int mode) const;
  };

  std::ostream &operator<<(std::ostream &str,const SF_Key &k);

  class Splitting_Function_Base;

  class SF_Lorentz {
  protected:

    ATOOLS::Flavour m_flavs[3], m_flspec;
    ATOOLS::Mass_Selector *&p_ms;

    SF_Coupling *p_cf;

    Splitting_Function_Base *p_sf;

    double m_zmin, m_zmax;
    int m_beam;

    static double s_pdfcut;

    double Lambda(const double &a,const double &b,const double &c);

  public:

    SF_Lorentz(const SF_Key &key);

    virtual ~SF_Lorentz();

    virtual double operator()
      (const double z,const double y,const double eta,
       const double scale,const double Q2,
       ATOOLS::Cluster_Amplitude *const sub=NULL) = 0;
    /**
      * Overestimate integrated over [zmin, zmax]
      */
    virtual double OverIntegrated(const double zmin,const double zmax,
				  const double scale,const double xbj) = 0;
    /**
      * Overestimate of operator
      */
    virtual double OverEstimated(const double z,const double y) = 0;
    /**
      * Generate z according to OverEstimated (solve OverIntegrated for new z)
      */
    virtual double Z() = 0;

    double JFF(const double &y) const;
    double JFI(const double &y,const double &eta,const double &scale,
	       const ATOOLS::Cluster_Amplitude *const sub) const;
    double JIF(const double &z,const double &y,
	       const double &eta,const double &scale,
	       const ATOOLS::Cluster_Amplitude *const sub) const;
    double JII(const double &z,const double &y,
	       const double &eta,const double &scale,
	       const ATOOLS::Cluster_Amplitude *const sub) const;

    inline const ATOOLS::Flavour &FlA() const { return m_flavs[0]; }
    inline const ATOOLS::Flavour &FlB() const { return m_flavs[1]; }
    inline const ATOOLS::Flavour &FlC() const { return m_flavs[2]; }

    inline const ATOOLS::Flavour &FlSpec() const { return m_flspec; }

    inline void SetFlSpec(const ATOOLS::Flavour &s) { m_flspec=s; }

    inline void SetBeam(const int beam) { m_beam=beam; }

    inline void SetSF(Splitting_Function_Base *const sf) { p_sf=sf; }

    inline const ATOOLS::Mass_Selector *MS() const { return p_ms; }

  };

  typedef ATOOLS::Getter_Function<SF_Lorentz,SF_Key,
				  std::less<std::string> > SFL_Getter;

  class SF_Coupling {
  protected:

    SF_Lorentz *p_lf;

    cstp::code m_type;
    double     m_cplfac;
    int        m_kfmode;

    /**
      * m_cplmax[i] = max_{scale} Coupling(scale, i)
      */
    std::vector<double> m_cplmax;

  public:

    SF_Coupling(const SF_Key &key);

    virtual ~SF_Coupling();

    virtual bool SetCoupling(MODEL::Model_Base *md,const double &k0sq,
			     const double &isfac,const double &fsfac) = 0;
    virtual double Coupling(const double &scale,const int pol,
			    ATOOLS::Cluster_Amplitude *const sub) = 0;
    virtual bool AllowSpec(const ATOOLS::Flavour &fl) = 0;

    virtual double CplFac(const double &scale) const;

    inline double MaxCoupling(const int mode) const { return m_cplmax[mode]; }

    inline void SetLF(SF_Lorentz *const lf) { p_lf=lf; }

  };

  class Shower;

  typedef ATOOLS::Getter_Function<SF_Coupling,SF_Key,
				  std::less<std::string> > SFC_Getter;

  struct SFC_Filler_Key {
    const MODEL::Model_Base *p_md;
    std::vector<SFC_Getter*> *p_gets;
    inline SFC_Filler_Key(const MODEL::Model_Base *const md,
			  std::vector<SFC_Getter*> *const gets):
      p_md(md), p_gets(gets) {}
  };

  typedef ATOOLS::Getter_Function
  <void,SFC_Filler_Key,std::less<std::string> > SFC_Filler_Getter;

  class Splitting_Function_Base {
  protected:
    SF_Lorentz  *p_lf;
    SF_Coupling *p_cf;
    cstp::code m_type;
    double m_lastint, m_symf, m_polfac, m_lpdf, m_rbmax;
    int m_on, m_bwon, m_qcd;
    std::vector<Parton*> m_specs;
    Parton *p_spec;
    PDF::PDF_Base **p_pdf;
    ATOOLS::RB_Key m_rbkey;
    double BWFactor(const double&,const double&,
		    const double&,const double&,
		    const double&) const;
  public:

    Splitting_Function_Base();
    Splitting_Function_Base(const SF_Key &key);

    virtual ~Splitting_Function_Base();

    virtual double operator()
      (const double z,const double y,const double eta,
       const double scale,const double Q2,
       ATOOLS::Cluster_Amplitude *const sub=NULL);
    virtual double OverIntegrated(const double,const double,
				  const double,const double);
    virtual double Overestimated(const double,const double);
    virtual double RejectionWeight(const double,const double,
				   const double,const double,const double);
    virtual double Z();         
    virtual double Phi(double z) const;

    virtual const ATOOLS::Flavour & GetFlavourA() const;
    virtual const ATOOLS::Flavour & GetFlavourB() const;
    virtual const ATOOLS::Flavour & GetFlavourC() const;
    virtual const ATOOLS::Flavour & GetFlavourSpec() const;

    virtual void ResetLastInt();
    virtual void ClearSpecs();

    virtual void SetRB(ATOOLS::RB_Map *const rbmap);

    double MEPSWeight(const double &z,const double &y,const double &eta,
		      const double &scale,const double &Q2,
		      const ATOOLS::Cluster_Amplitude *const sub) const;

    double GetXPDF(const double &scale,const double &x,
		   const ATOOLS::Flavour &a,const int beam,
		   const int mode=0);

    Parton *SelectSpec();
    bool PureQCD() const;

    Parton *SetSpec(Parton *const spec);

    inline double     Last() const    { return m_lastint; } 
    inline cstp::code GetType() const { return m_type;    } 

    inline void SetFlavourSpec(ATOOLS::Flavour flav) { p_lf->SetFlSpec(flav); }
    inline void AddSpec(Parton *const s) { m_specs.push_back(s); }

    inline Parton *Spec() const { return p_spec; }

    inline int On() const { return m_on; }
    template <class _Type> inline _Type *Get() { return dynamic_cast<_Type*>(this); }

    inline SF_Lorentz  *Lorentz() const  { return p_lf; }
    inline SF_Coupling *Coupling() const { return p_cf; }

    inline void SetPDF(PDF::PDF_Base **const pdf) { p_pdf=pdf; }

    inline double RBMax() const { return m_rbmax; }
    inline double SymFac() const { return m_symf; }

    inline void SetRB(const double &rbmax) { m_rbmax=rbmax; }

    friend std::ostream& operator<<(std::ostream &,const Splitting_Function_Base &);
  };
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
