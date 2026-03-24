#ifndef METOOLS_Current_Library_FF_0_PP_H
#define METOOLS_Current_Library_FF_0_PP_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  enum FF_0_PP_mode {
    pipi_plus       = 1,
    KK_plus         = 2,
    etapi_plus      = 3,
    etaprimepi_plus = 4,
    Kpi_plus        = 11,
    Keta_plus       = 12,
    Ketaprime_plus  = 13,
    pipi_0          = 101,
    KK_0            = 102,
    etapi_0         = 103,
    etaprimepi_0    = 104,
    unknown         = 999
  };
  
  class FF_0_PP_Base : public FormFactor_Base {
  protected:
    FF_0_PP_mode        m_mode;
    double              m_norm;
    Summed_Propagator * p_props;

    void    FixMode();
    virtual Complex FF_RChiPT(const double & Q2);
  public:
    FF_0_PP_Base(const FF_Parameters & params) :
      FormFactor_Base(params),
      m_mode(FF_0_PP_mode::unknown),
      m_norm(0.), p_props(NULL) {
      m_ffmodel = ff_model::RChiPT;
      FixMode();
    }
    ~FF_0_PP_Base() { if (p_props) delete p_props; }

    Complex operator()(const ATOOLS::Vec4D_Vector& moms);
    inline const double & GetNorm() const { return m_norm; }
  };
}

#endif
