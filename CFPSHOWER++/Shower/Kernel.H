#ifndef CFPSHOWER_Shower_Kernel_H
#define CFPSHOWER_Shower_Kernel_H

#include "CFPSHOWER++/Shower/Splitting.H"
#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Calculators/Gauge_Base.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Math/MathTools.H"
#include <iostream>

namespace PDF {
  class PDF_Base;
}

namespace CFPSHOWER {
  class Kernel {
  private:
    kernel_type::code m_type;
    bool              m_swapped;
    ATOOLS::Flavour_Vector        m_flavs;
    const ATOOLS::Mass_Selector * p_msel;
    
    SF_Base       * p_sf; 
    Gauge_Base    * p_gauge;
    PDF::PDF_Base * p_pdf[2];
    double          m_pdfminvalue, m_pdfxmin;

    Weight * MakeWeight(const Splitting & split,const double & overfac=1.);
  public:
    Kernel(const Kernel_Info & info);
    inline ~Kernel() {
      delete p_sf;
      delete p_gauge;
    }
    
    double Integral(Splitting & split,const ATOOLS::Mass_Selector * msel);
    bool   Generate(Splitting & split,const ATOOLS::Mass_Selector * msel,
		    const double & overfac=1.);
    bool   FillOffsprings(Splitting & split);

    double GetXPDF(const double & x,const double & Q2,
		   const ATOOLS::Flavour & flav,const size_t beam);
    
    inline void SetSF(SF_Base * sf)                            { p_sf    = sf; }
    inline void SetGauge(Gauge_Base * gaugepart)               { p_gauge = gaugepart; }
    inline void SetMSel(const ATOOLS::Mass_Selector * ms)      { p_msel  = ms; }
    inline void SetPDF(const size_t beam, PDF::PDF_Base * pdf) { p_pdf[beam]   = pdf; }
    inline void SetPDFMinValue(const double & value)           { m_pdfminvalue = value; }
    inline void SetPDFXMin(const double & xmin)                { m_pdfxmin     = xmin; }
    
    inline SF_Base                      * GetSF()    const { return p_sf; } 
    inline Gauge_Base                   * GetGauge() const { return p_gauge; }
    inline const ATOOLS::Flavour_Vector & GetFlavs() const { return m_flavs; }
    inline const ATOOLS::Mass_Selector  * GetMSel()  const { return p_msel; }
    inline PDF::PDF_Base * GetPDF(const size_t beam) const { return p_pdf[beam]; }
    inline const double  & PDFMinValue()             const { return m_pdfminvalue; }
    inline const double  & PDFXMin()                 const { return m_pdfxmin; }
  };
  
  std::ostream &operator<<(std::ostream &s,const Kernel & kernel);

  typedef ATOOLS::Getter_Function<Kernel,Kernel_Info,std::less<std::string> > Kernel_Getter;
  typedef std::vector<Kernel *> Kernel_Vector;
}

#endif
