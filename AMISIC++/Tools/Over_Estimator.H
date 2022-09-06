#ifndef AMISIC_Tools_Over_Estimator_H
#define AMISIC_Tools_Over_Estimator_H

#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"

namespace AMISIC {
  class MI_Processes;
  
  class Over_Estimator {
  private:
    size_t m_nbins;
    double m_muR_fac, m_muF_fac;
    double m_pref, m_bfac, m_xsnd, m_pt02, m_ptmin2, m_s, m_xt;
    PDF::PDF_Base         * p_pdf[2];
    double                  m_xmin[2];
    MODEL::Running_AlphaS * p_alphaS;
    MI_Processes *          p_procs;

    void   FixMaximum();
    double ApproxME(const double & pt2);
    double ExactME(const double & pt2);
  public:
    Over_Estimator();
    ~Over_Estimator();

    void   Initialize(MI_Processes * procs);
    void   Update();

    double operator()(const double & pt2);
    double TrialPT2(const double & Q2);
    
    inline void SetPref(const double & pref) { m_pref = pref; }
    inline void SetBFac(const double & bfac) { m_bfac = bfac; }
    inline void SetXSnd(const double & xsnd) { m_xsnd = xsnd; }
    inline void SetPT02(const double & pt02) { m_pt02 = pt02; }

    void Test(const double & Q2,const long int & n);
  };
}

#endif
