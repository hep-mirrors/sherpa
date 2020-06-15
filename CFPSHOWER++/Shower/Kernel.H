#ifndef CFPSHOWER_Shower_Kernel_H
#define CFPSHOWER_Shower_Kernel_H

#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Calculators/Gauge_Base.H"
#include "CFPSHOWER++/Calculators/Kinematics_Base.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "CFPSHOWER++/Tools/Configuration.H"
#include "CFPSHOWER++/Tools/Splitting.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Math/MathTools.H"
#include <iostream>

namespace PDF {
  class PDF_Base;
}

namespace CFPSHOWER {
  class Kernel {
  private:
    const ATOOLS::Mass_Selector * p_msel;
    
    SF_Base         * p_sf; 
    Gauge_Base      * p_gauge;
    Kinematics_Base * p_kinematics;
    PDF::PDF_Base   * p_pdf[2];
    double            m_pdfminvalue, m_pdfxmin, m_enhance;

    Weight * MakeWeight(const Splitting & split,const double & overfac=1.);
  public:
    Kernel(const Kernel_Info & info);
    inline ~Kernel() {
      delete p_sf;
      delete p_gauge;
    }
    
    double Integral(Splitting & split,const ATOOLS::Mass_Selector * msel);
    bool   Generate(Splitting & split,Configuration & config,
		    const ATOOLS::Mass_Selector * msel,const double & overfac=1.);
    bool   UpdateSystem(Splitting & split,Configuration & config);
    bool   FillOffsprings(Splitting & split);

    double GetXPDF(const double & x,const double & Q2,
		   const ATOOLS::Flavour & flav,const size_t beam);
    
    inline void SetSF(SF_Base * sf)                            { p_sf          = sf; }
    inline void SetGauge(Gauge_Base * gaugepart)               { p_gauge       = gaugepart; }
    inline void SetKinematics(Kinematics_Base * kinpart)       { p_kinematics  = kinpart; }
    inline void SetMSel(const ATOOLS::Mass_Selector * ms)      { p_msel        = ms; }
    inline void SetPDF(const size_t beam, PDF::PDF_Base * pdf) { p_pdf[beam]   = pdf; }
    inline void SetPDFMinValue(const double & value)           { m_pdfminvalue = value; }
    inline void SetPDFXMin(const double & xmin)                { m_pdfxmin     = xmin; }
    inline void SetEnhanceFactor(const double & enhance)       { m_enhance     = enhance; }

    inline SF_Base                      * GetSF()         const { return p_sf; }
    inline Kinematics_Base              * GetKinematics() const { return p_kinematics; }
    inline const kernel_type::code      & Type()          const { return p_sf->Type(); }
    inline const log_type::code         & LogType()       const { return p_sf->LogType(); }
    inline Gauge_Base                   * GetGauge()      const { return p_gauge; }
    inline const ATOOLS::Flavour        & GetSplit()      const { return p_sf->GetSplit(); }
    inline const ATOOLS::Flavour_Vector & GetFlavs()      const { return p_sf->GetFlavs(); }
    inline const std::vector<size_t>    & Tags()          const { return p_sf->Tags(); }
    inline const size_t & Tags(const size_t & pos)        const { return p_sf->Tags(pos); }
    inline const ATOOLS::Mass_Selector  * GetMSel()       const { return p_msel; }
    inline PDF::PDF_Base * GetPDF(const size_t beam)      const { return p_pdf[beam]; }
    inline const double & PDFMinValue()                   const { return m_pdfminvalue; }
    inline const double & PDFXMin()                       const { return m_pdfxmin; }
    inline const double & EnhanceFactor()                 const { return m_enhance; }  
  };
  
  std::ostream &operator<<(std::ostream &s,const Kernel & kernel);

  typedef ATOOLS::Getter_Function<Kernel,Kernel_Info,std::less<std::string> > Kernel_Getter;
  typedef std::vector<Kernel *> Kernel_Vector;
}

#endif
