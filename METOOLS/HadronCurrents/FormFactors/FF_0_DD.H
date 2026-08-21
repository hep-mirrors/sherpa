#ifndef METOOLS_Current_Library_FF_0_DD_H
#define METOOLS_Current_Library_FF_0_DD_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  enum FF_0_DD_mode {
    pp              = 1,
    nn              = 2,
    unknown         = 999
  };
  
  class FF_0_DD_Base : public FormFactor_Base {
  protected:
    FF_0_DD_mode m_mode;
    double       m_norm;

    void         FixMode();
  public:
    FF_0_DD_Base(const FF_Parameters & params) :
      FormFactor_Base(params),
      m_mode(FF_0_DD_mode::unknown), m_norm(0.)
    { FixMode(); }
    ~FF_0_DD_Base() {}

    virtual Complex operator()(const ATOOLS::Vec4D_Vector& moms) = 0;
    inline const double & GetNorm() const { return m_norm; }
  };
}

#endif
