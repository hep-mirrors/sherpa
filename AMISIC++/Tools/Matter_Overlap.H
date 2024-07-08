#ifndef AMISIC_Tools_Matter_Overlap_H
#define AMISIC_Tools_Matter_Overlap_H

#include "AMISIC++/Tools/Lookup_Tables.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Vector.H"


namespace AMISIC {
  class MO_Integrand;

  class Matter_Overlap : public ATOOLS::Function_Base {
  private:
    bool m_dynamic;
    REMNANTS::Form_Factor * p_ffs[2];
    REMNANTS::matter_form   m_form[2];
    double m_bstep, m_bmax;
    axis * p_bbins;
    double m_integral, m_kradius, m_norm, m_invGeV2fm;
    double m_radius[4], m_radius2[4], m_rnorm[4], m_fraction[4];
double m_fixradius, m_dynradius2, m_maxradius;
    double m_xmin[2], m_xmax[2];


    void   InitializeFFParams(PDF::ISR_Handler * const isr);
    void   InitializeStaticFFParams();
    void   InitializeDynamicFFParams(PDF::ISR_Handler * const isr);
    double CalculateIntegral();
    void   Output(const double & check=0.);
 public:
    Matter_Overlap();
    ~Matter_Overlap();

    void   Initialize(REMNANTS::Remnant_Handler * const rh,
		      PDF::ISR_Handler * const isr);
    void   FixDynamicRadius(const double & x1=-1.,  const double & x2=-1.,
			    const double & Q21=-1., const double & Q22=-1.);
    double operator()(double b);
    double MaxValue(const double & b);

    double SelectB() const;

    ATOOLS::Vec4D
    SelectPositionForScatter(const double & b,
			     const double & x0, const double & Q20,
			     const double & x1, const double & Q21) const;

    inline void SetKRadius(const double & k)   { m_kradius = k; }
    inline void SetBmax(const double & bmax)   { m_bmax    = bmax; }
    inline void SetBAxis(axis * bbins)         { p_bbins   = bbins; }
    inline const bool   & IsDynamic()      const { return m_dynamic; }
    inline const double & Bmax()           const { return m_bmax; }
    inline const double & Bstep()          const { return m_bstep; }
    inline const double & Integral()       const { return m_integral; }
    inline const double & KRadius()        const { return m_kradius; }
    inline const double & DynamicRadius2() const { return m_dynradius2; }
    inline const double & FixRadius()      const { return m_fixradius; }
    inline       axis   * GetBBins()       const { return p_bbins; }
  };

  class MO_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
  public:
    MO_Integrand(Matter_Overlap * mo) : p_mo(mo) {}
    ~MO_Integrand() = default;
    double operator()(double b);
  };
}

#endif
