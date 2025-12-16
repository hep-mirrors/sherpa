#ifndef REMNANTS_Tools_Pseudo_TMD_H
#define REMNANTS_Tools_Pseudo_TMD_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace REMNANTS {
  class Pseudo_TMD {
  private:
    PDF::PDF_Base *            p_pdf;
    std::list<ATOOLS::Flavour> m_pdfpartons;
    ATOOLS::Flavour            m_bunch;
    double m_xmin,m_xmax,m_pt2min,m_pt2max;
    double m_Q02,m_Y0,m_sigma02,m_x0;                         
    double m_x,m_Q2,m_pt2,m_Y;
    std::map<ATOOLS::Flavour,double> m_xpdfmax;
    
    const double Sigma2(const double & y) const;
    void Test();
  public:
    Pseudo_TMD(PDF::PDF_Base * pdf,const double & Y0,const double & sigma0);
    ~Pseudo_TMD();

    void   Scan();
    void   Calculate(const double & x,const double & Q2,
		     const double & pt2,const double & Y);
    double AllPartons(const double & x,const double & Q2,
		      const double & pt2,const double & Y);
    double XPDF(const ATOOLS::Flavour & flav,const bool & defmax=false);
    double XTMD(const ATOOLS::Flavour & flav,const bool & defmax=false);

    inline const std::list<ATOOLS::Flavour> & GetFlavours() const {
      return m_pdfpartons;
    } 
    inline const double &          XMin()  const { return m_xmin; }
    inline const double &          XMax()  const { return m_xmax; }
    inline double                  Q2Min() const { return 0.; }
    inline double                  YMin()  const { return 0.; }
    inline const ATOOLS::Flavour & Bunch() const { return m_bunch; }
  };

  class Pseudo_TMD_PT2Test : public ATOOLS::Function_Base {
    Pseudo_TMD * p_pdf;
    double       m_x,m_Q2,m_y;
    ATOOLS::Flavour m_flav;
  public:
    Pseudo_TMD_PT2Test(Pseudo_TMD * pdf) :
      p_pdf(pdf), m_Q2(0.) { }
    ~Pseudo_TMD_PT2Test() {}
    double operator()(double pt2);
    inline void SetXQ2(const double & x,const double & Q2) { m_x=x;m_Q2=Q2; }
    inline void Sety(const double & y)                     { m_y=y; }
    inline void SetFlav(const ATOOLS::Flavour & flav=
			ATOOLS::Flavour(kf_none))          { m_flav=flav; }  
  };
}

#endif
