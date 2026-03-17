#ifndef METOOLS_Current_Library_FF_0_PP_H
#define METOOLS_Current_Library_FF_0_PP_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  class FF_0_PP_Base : public FormFactor_Base {
  protected:
    enum class PSmode {
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
      Kpi_0           = 111,
      Keta_0          = 112,
      Ketaprime_0     = 113,
      unknown         = 999
    };
  protected:
    PSmode m_PSmode;
    double m_norm;

    Summed_Propagator * p_props;
    
    void FixMode();
    void FixNorm();
    virtual void Construct() {};
  public:
    FF_0_PP_Base(const FF_Parameters & params);
    ~FF_0_PP_Base();
    Complex operator()(const ATOOLS::Vec4D_Vector& moms);
    inline const double & GetNorm() const { return m_norm; }
  };
  
  class Fplus_0_PP : public FF_0_PP_Base {
    void Construct() override;
  public:
    Fplus_0_PP(const FF_Parameters & params);
    ~Fplus_0_PP() {}
  };

  class F0_0_PP : public FF_0_PP_Base {
    void Construct() {};
  public:
    F0_0_PP(const FF_Parameters & params);
    ~F0_0_PP() {}
  };
}

#endif
