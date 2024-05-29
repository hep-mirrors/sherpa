#ifndef AMISIC_Tools_Over_Estimator_H
#define AMISIC_Tools_Over_Estimator_H

#include "AMISIC++/Tools/Lookup_Tables.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include <array>

namespace PDF {
  class ISR_Handler;
  class PDF_Base;
}

namespace AMISIC {
  class MI_Processes;

  class Over_Estimator {
  private:
    double m_muR_fac, m_muF_fac, m_pt02, m_ptmin2;
    std::array<double, 2> m_xmin, m_xmax;
    double m_s, m_pref, m_bfac, m_xsnd;
    size_t m_pt2bins;

    std::array<PDF::PDF_Base *, 2> p_pdf;
    MODEL::One_Running_AlphaS * p_alphaS;
    OneDim_Table              * p_prefs;

    void   FixMaximum(MI_Processes * procs,axis * sbins,axis * pt2bins);
    double ApproxME(const double & pt2,const double & xt);
    double ExactME(const double & pt2,const double & xt);

    void Output();
  public:
    Over_Estimator();
    ~Over_Estimator();

    void   Initialize(PDF::ISR_Handler * isr,MI_Processes * procs,
		      axis * sbins,axis * pt2bins);
    void   UpdateS(const double & s,const double & xsnd, const double & pt02);

    double operator()(const double & pt2, const double & yvol);
    double TrialPT2(const double & Q2);

    inline const double & BFac() const { return m_bfac; }
    
    inline void SetBFac(const double & bfac) { m_bfac = bfac; }
    inline void SetXSnd(const double & xsnd) { m_xsnd = xsnd; }

    void Test(const double & Q2,const long int & n);
  };
}

#endif
