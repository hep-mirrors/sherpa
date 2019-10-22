#ifndef MCATNLO_Showers_SF_Lorentz_H
#define MCATNLO_Showers_SF_Lorentz_H

#include "MCATNLO/Showers/SF_Key.H"

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Flavour.H"
#include "EXTAMP/CS_Dipole.H"

namespace ATOOLS {
  class Cluster_Amplitude;
}

namespace MCATNLO {

  class SF_Coupling;
  class Splitting_Function_Base;

  class SF_Lorentz {
  protected:

    ATOOLS::Flavour m_flavs[3], m_flspec;
    const ATOOLS::Mass_Selector *&p_ms;

    SF_Coupling *p_cf;

    Splitting_Function_Base *p_sf;

    double m_zmin, m_zmax;
    int m_beam, m_col;

    double m_lastJ;

    std::pair<double, double> m_pdfmin;

    static double s_kappa;

    double Lambda(const double &a,const double &b,const double &c) const;
    bool PDFValueAllowedAsDenominator(const double& val,
                                      const double& eta);

  public:

    SF_Lorentz(const SF_Key &key);

    virtual ~SF_Lorentz();

    virtual double Scale(const double z,const double y,
			 const double scale,const double Q2) const;

    virtual double operator()
      (const double z,const double y,const double eta,
       const double scale,const double Q2,
       ATOOLS::Cluster_Amplitude *const sub=NULL) = 0;
    virtual double AsymmetryFactor(const double z,const double y,
				   const double Q2);
    /**
      * Overestimate integrated over [zmin, zmax]
      */
    virtual double OverIntegrated(const double zmin,const double zmax,
				  const double scale,const double xbj) = 0;
    /**
      * Overestimate of operator
      */
    virtual double OverEstimated(const double z,const double y,const double phi=0) = 0;
    /**
      * Generate z according to OverEstimated (solve OverIntegrated for new z)
      */
    virtual double Z() = 0;

    double LastJ() const { return m_lastJ; }
    void SetLastJ(double J) { m_lastJ = J; }

    double JFF(const double &y,const double &mui2,const double &muj2,
	       const double &muk2,const double &muij2);
    double JFI(const double &y,const double &eta,const double &scale,
	       const ATOOLS::Cluster_Amplitude *const sub);
    double JIF(const double &z,const double &y,
	       const double &eta,const double &scale,
	       const ATOOLS::Cluster_Amplitude *const sub);
    double JII(const double &z,const double &y,
	       const double &eta,const double &scale,
	       const ATOOLS::Cluster_Amplitude *const sub);

    inline const ATOOLS::Flavour &FlA() const { return m_flavs[0]; }
    inline const ATOOLS::Flavour &FlB() const { return m_flavs[1]; }
    inline const ATOOLS::Flavour &FlC() const { return m_flavs[2]; }

    inline const ATOOLS::Flavour &FlSpec() const { return m_flspec; }

    inline void SetFlSpec(const ATOOLS::Flavour &s) { m_flspec=s; }

    inline int GetBeam() const { return m_beam; }
    inline void SetBeam(const int beam) { m_beam=beam; }

    inline void SetSF(Splitting_Function_Base *const sf) { p_sf=sf; }

    inline const ATOOLS::Mass_Selector *MS() const { return p_ms; }

    inline int Col() const { return m_col; }

    inline static void SetKappa(const double &kap) { s_kappa=kap; }

    EXTAMP::DipoleCase m_dipole_case;

  };

  typedef ATOOLS::Getter_Function<SF_Lorentz,SF_Key,
				  std::less<std::string> > SFL_Getter;

}

#endif
