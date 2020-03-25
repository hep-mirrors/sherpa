#ifndef MCATNLO_Showers_Splitting_Function_Base_H
#define MCATNLO_Showers_Splitting_Function_Base_H

#include "MCATNLO/Showers/SF_Key.H"
#include "MCATNLO/Showers/SF_Coupling.H"
#include "MCATNLO/Showers/SF_Lorentz.H"

namespace ATOOLS { 
  class Cluster_Amplitude;
}

namespace PDF { class PDF_Base; }

namespace MCATNLO {

  class Parton;
  class Shower;

  class Splitting_Function_Base {
  private:
    void RegisterDefaults() const;
  protected:
    SF_Lorentz  *p_lf;
    SF_Coupling *p_cf;
    cstp::code m_type;
    double m_lastint, m_lastacceptwgt, m_lastscale;
    double m_symf, m_polfac, m_lpdf, m_mth;
    int m_on, m_qcd, m_scs;
    std::vector<Parton*> m_specs;
    std::vector<double> m_lastints;
    Parton *p_spec;
    PDF::PDF_Base **p_pdf;
  public:

    Splitting_Function_Base();
    Splitting_Function_Base(const SF_Key &key);

    virtual ~Splitting_Function_Base();

    virtual double operator()
      (const double z,const double y,const double eta,
       const double scale,const double Q2,
       const Color_Info &ci=Color_Info(),
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

    virtual int GetCol() const;

    virtual void ResetLastInt();
    virtual void ClearSpecs();

    double AsymmetryFactor(const double z,const double y,
			   const double Q2);

    double MEPSWeight(const double &z,const double &y,const double &eta,
		      const double &scale,const double &Q2,
		      const ATOOLS::Cluster_Amplitude *const sub) const;

    double GetXPDF(const double &scale,const double &x,
		   const ATOOLS::Flavour &a,const int beam,
		   const int mode=0);

    Parton *SelectSpec();
    bool PureQCD() const;

    Parton *SetSpec(Parton *const spec);

    void ColorPoint(Parton *const p) const;

    double ColorWeight(Parton *const p) const;

    inline double     Last()                 const { return m_lastint;       } 
    virtual double    LastAcceptanceWeight() const { return m_lastacceptwgt; }
    virtual double    LastScale()            const { return m_lastscale;     }
    inline cstp::code GetType()              const { return m_type;          } 

    inline void SetFlavourSpec(ATOOLS::Flavour flav) { p_lf->SetFlSpec(flav); }
    inline void AddSpec(Parton *const s) { m_specs.push_back(s); }

    inline Parton *Spec() const { return p_spec; }

    inline int On() const { return m_on; }
    template <class _Type> inline _Type *Get() { return dynamic_cast<_Type*>(this); }

    inline SF_Lorentz  *Lorentz() const  { return p_lf; }
    inline SF_Coupling *Coupling() const { return p_cf; }

    PDF::PDF_Base** PDF() { return p_pdf; }
    inline void SetPDF(PDF::PDF_Base **const pdf) { p_pdf=pdf; }

    inline double SymFac() const { return m_symf; }

    inline void SetMassThreshold(const double &mth) { m_mth=mth; }

    inline void SetScaleScheme(const int &scs) { m_scs=scs; }

    inline int ScaleScheme() const { return m_scs; }

    friend std::ostream& operator<<(std::ostream &,const Splitting_Function_Base &);
  };
}

#endif
