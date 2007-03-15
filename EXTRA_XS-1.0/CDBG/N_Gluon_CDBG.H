#ifndef N_Gluon_CDBG_H
#define N_Gluon_CDBG_H

#include "CDBG_Amplitude.H"

namespace PHASIC {

  class Color_Integrator;
  class Helicity_Integrator;

}

namespace EXTRAXS {

  class N_Gluon_CDBG {
  private:

    size_t m_nin, m_nout;

    CDBG_Amplitude m_ampl;

    PHASIC::Color_Integrator    *p_colint;
    PHASIC::Helicity_Integrator *p_helint;

    std::vector<ATOOLS::Vec4D> m_moms;
  
    double m_as, m_sf;

    int m_mode, m_tests;

    double Factorial(const double &n);

  public:

    // constructor
    N_Gluon_CDBG(const size_t nin,const size_t nout,
		 const std::vector<ATOOLS::Flavour> &flavs);

    // destructor
    ~N_Gluon_CDBG();

    // member functions
    double Differential(const std::vector<ATOOLS::Vec4D> &momenta);

    bool GaugeTest(std::vector<ATOOLS::Vec4D> momenta);

    // inline functions
    inline void SetMode(const int &mode)    { m_mode=mode;   }
    inline void SetTests(const int &tests)  { m_tests=tests; }

    inline void SetColorIntegrator(PHASIC::Color_Integrator *const cint)
    { p_colint=cint; }
    inline void SetHelicityIntegrator(PHASIC::Helicity_Integrator *const hint)
    { p_helint=hint; }

  };// end of class N_Gluon_CDBG

}// end of namespace EXTRAXS

#endif

