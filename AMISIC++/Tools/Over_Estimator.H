#ifndef AMISIC_Tools_Over_Estimator_H
#define AMISIC_Tools_Over_Estimator_H

#include "AMISIC++/Tools/Lookup_Tables.H"
#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"

namespace AMISIC {
  class MI_Processes;

  class Over_Estimator {
  private:
    double m_muR_fac, m_muF_fac, m_pt02, m_ptmin2;
    std::array<double, 2> m_xmin;
    double m_s, m_pref, m_bfac, m_xsnd;
    size_t m_pt2bins;

    std::array<PDF::PDF_Base *, 2> p_pdf;
    MODEL::Running_AlphaS * p_alphaS;
    MI_Processes          * p_procs;
    OneDim_Table          * p_prefs;

    void   FixMaximum();
    double ApproxME(const double & pt2,const double & xt);
    double ExactME(const double & pt2,const double & xt);
  public:
    Over_Estimator();
    ~Over_Estimator();

    void   Initialize(MI_Processes * procs);
    void   UpdateS();

    double operator()(const double & pt2, const double & yvol);
    double TrialPT2(const double & Q2);

    inline void SetBFac(const double & bfac) { m_bfac = bfac; }
    inline void SetXSnd(const double & xsnd) { m_xsnd = xsnd; }

    void Test(const double & Q2,const long int & n);
  };
}

#endif
