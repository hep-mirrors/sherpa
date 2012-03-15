#ifndef SHRIMPS_Beam_Remnants_Continued_PDF_H
#define SHRIMPS_Beam_Remnants_Continued_PDF_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SHRIMPS {
  class PDF_Kernel_Base : public ATOOLS::Function_Base {
  protected:
    PDF::PDF_Base              * p_pdf;
    ATOOLS::Flavour              m_bunch;
    std::list<ATOOLS::Flavour> * p_pdfpartons;
    double                       m_xmin,m_xmax,m_Q02;
  public:
    PDF_Kernel_Base(PDF::PDF_Base * pdf,
		    const ATOOLS::Flavour & bunch,
		    std::list<ATOOLS::Flavour> * pdfpartons,
		    const double & Q02=-1.) :
      p_pdf(pdf), m_bunch(bunch), p_pdfpartons(pdfpartons), 
      m_xmin(p_pdf->XMin()), m_xmax(p_pdf->XMax()), 
      m_Q02(ATOOLS::Max(Q02,p_pdf->Q2Min())) { }
    virtual ~PDF_Kernel_Base() {}
    virtual double operator()(double x) = 0;
  };
  
  class Continued_PDF {
    class Sea_Kernel : public PDF_Kernel_Base {
    public:
      Sea_Kernel(PDF::PDF_Base * pdf,
		 const ATOOLS::Flavour & bunch,
		 std::list<ATOOLS::Flavour> * pdfpartons,
		 const double & Q02=-1.) :
	PDF_Kernel_Base(pdf,bunch,pdfpartons,Q02) {}
      double operator()(double x);
    };
    class Valence_Kernel : public PDF_Kernel_Base {
    public:
      Valence_Kernel(PDF::PDF_Base * pdf,
		     const ATOOLS::Flavour & bunch,
		     std::list<ATOOLS::Flavour> * pdfpartons,
		     const double & Q02=-1.) :
	PDF_Kernel_Base(pdf,bunch,pdfpartons,Q02) {}
      double operator()(double x);
    };

  private:
    PDF::PDF_Base                  * p_pdf;
    std::list<ATOOLS::Flavour>       m_pdfpartons;
    ATOOLS::Flavour                  m_bunch;
    double                           m_xmin,m_xmax,m_Q02,m_Vnorm,m_Snorm;
    double                           m_x,m_Q2;
    std::map<ATOOLS::Flavour,double> m_xpdfmax;
    std::map<ATOOLS::Flavour,double> m_xmaxpdf;

    void   CalculateNorms();
  public:
    Continued_PDF(PDF::PDF_Base * pdf,const ATOOLS::Flavour & bunch);
    ~Continued_PDF();

    const double XPDF(const ATOOLS::Flavour & flav,const bool & defmax=false);

    inline void Calculate(const double & x,const double & Q2) {
      m_x  = x; 
      m_Q2 = Q2;
      if (Q2<m_Q02) p_pdf->Calculate(m_x,m_Q02);
      else p_pdf->Calculate(m_x,m_Q2);
    }
    inline double AllPartons(const double & x,const double & Q2) {
      Calculate(x,Q2);
      double val(0.), test(0.);
      for (std::list<ATOOLS::Flavour>::iterator flit=m_pdfpartons.begin();
	   flit!=m_pdfpartons.end();flit++) {
	val += test = XPDF((*flit));
	if (test>m_xpdfmax[(*flit)]) m_xpdfmax[(*flit)] = test;
      }
      return val;
    }
    inline const double &          XMin()  const { return m_xmin; }
    inline const double &          XMax()  const { return m_xmax; }
    inline const double            Q2Min() const { return 0.; }
    inline const ATOOLS::Flavour & Bunch() const { return m_bunch; }

    inline const double   XPDFMax(const ATOOLS::Flavour & flav) {
      ATOOLS::Flavour flav1(flav);
      if (flav.IsDiQuark()) {
	flav1 = ATOOLS::Flavour(kf_u);
	if (m_bunch==ATOOLS::Flavour(kf_p_plus).Bar() && flav.IsAnti()) 
	  flav1 = flav1.Bar();
      }
      std::map<ATOOLS::Flavour,double>::iterator found(m_xpdfmax.find(flav1));
      if (found==m_xpdfmax.end()) return -1.;
      return found->second;
    }
    inline const double   XMaxPDF(const ATOOLS::Flavour & flav) {
      ATOOLS::Flavour flav1(flav);
      if (flav.IsDiQuark()) {
	flav1 = ATOOLS::Flavour(kf_u);
	if (m_bunch==ATOOLS::Flavour(kf_p_plus).Bar() && flav.IsAnti()) 
	  flav1 = flav1.Bar();
      }
      std::map<ATOOLS::Flavour,double>::iterator found(m_xmaxpdf.find(flav1));
      if (found==m_xmaxpdf.end()) return -1.;
      return found->second;
    }
  };
}

#endif
