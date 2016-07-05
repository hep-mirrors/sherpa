#ifndef PDF_Main_NLOMC_Base_H
#define PDF_Main_NLOMC_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Flavour.H"

#include <string>

namespace SHERPA { class Variation_Weights; }

namespace ATOOLS {

  class Cluster_Amplitude;
  class Data_Reader;
  struct RB_Data;

}

namespace MODEL { class Model_Base; }

namespace PHASIC { class Process_Base; }

namespace PDF {

  class ISR_Handler;
  class Shower_Base;

  typedef std::pair<double,double> SH_Pair;

  struct Dipole_Value {
    PHASIC::Process_Base *p_proc;
    double m_d, m_b, m_t, m_mu2;
    bool m_trig;
    inline Dipole_Value(PHASIC::Process_Base *const proc=NULL,
			const double &d=0.0,const double &b=0.0,
			const double &t=0.0,const double &mu2=0.0,
			const bool trig=true):
      p_proc(proc), m_d(d), m_b(b), m_t(t), m_mu2(mu2), m_trig(trig) {}
  };// end of struct Dipole_Value

  class NLOMC_Base {
  protected:

    std::string m_name;

    Shower_Base *p_shower;

    double m_kt2min, m_weight;

    SHERPA::Variation_Weights * p_variationweights;

  public:

    NLOMC_Base(const std::string &name);

    virtual ~NLOMC_Base();

    virtual int GeneratePoint(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    virtual void AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    static void ShowSyntax(const int mode);

    inline const std::string &Name() const { return m_name; }

    inline double KT2Min() const { return m_kt2min; }
    inline double Weight() const { return m_weight; }

    inline void SetShower(Shower_Base *const shower) { p_shower=shower; }

    inline void SetVariationWeights(SHERPA::Variation_Weights * varwgts)
    { p_variationweights = varwgts; }

  };// end of class NLOMC_Base

  struct NLOMC_Key {
    MODEL::Model_Base *p_model;
    ISR_Handler *p_isr;
    ATOOLS::Data_Reader *p_read;
    inline NLOMC_Key(MODEL::Model_Base *const model,
		      ISR_Handler *const isr,
		      ATOOLS::Data_Reader *const read):
      p_model(model), p_isr(isr), p_read(read) {}
  };//end of struct NLOMC_Key

  typedef ATOOLS::Getter_Function
  <NLOMC_Base,NLOMC_Key> NLOMC_Getter;

}// end of namespace PDF

#endif
