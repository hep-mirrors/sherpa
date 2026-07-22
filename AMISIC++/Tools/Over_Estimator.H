#ifndef AMISIC_Tools_Over_Estimator_H
#define AMISIC_Tools_Over_Estimator_H

#include "ATOOLS/Math/Lookup_Tables.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include <array>

namespace PDF {
  class ISR_Handler;
  class PDF_Base;
}

namespace AMISIC {
  using ATOOLS::axis;
  using ATOOLS::OneDim_Table;

  class Interaction_Probability;
  class Matter_Overlap;
  class MI_Processes;
  class Matter_Overlap;

  class Over_Estimator {
  private:
    double m_muR_fac, m_muF_fac, m_pt02, m_ptmin2;
    std::array<double, 2> m_xmin, m_xmax;
    double m_s, m_pref;
    size_t m_n_variations;
    size_t m_npt2bins;

    std::array<PDF::PDF_Base *, 2> p_pdf;
    MODEL::One_Running_AlphaS * p_alphaS;
    OneDim_Table              * p_prefs;

    void   FixMaximum(MI_Processes * procs,axis * sbins,
                      Interaction_Probability * pint,Matter_Overlap * mo);
    double ApproxME(const double & pt2,const double & xt);

    void Output();
  public:
    Over_Estimator();
    ~Over_Estimator();

    void   Initialize(PDF::ISR_Handler * isr,MI_Processes * procs, axis * sbins,
                      Interaction_Probability * pint,Matter_Overlap * mo);
    void   UpdateS(const double & s,const double & pt02,const double & ptmin2);

    double operator()(const double & pt2,const double & yvol);
    double TrialPT2(const double & Q2);

    void SetVariations(size_t n_variations) { m_n_variations = n_variations; }
  };
}

#endif
