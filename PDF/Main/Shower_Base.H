#ifndef PDF_Main_Shower_Base_H
#define PDF_Main_Shower_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/Event_Weights.H"
#include "ATOOLS/Phys/Flavour.H"

#include <string>

namespace ATOOLS {
  class Cluster_Amplitude;
  class Blob;
}

namespace MODEL    { class Model_Base; }

namespace REMNANTS { class Remnant_Handler; }

namespace PDF {

  class ISR_Handler;
  class Cluster_Definitions_Base;

  class Shower_Base {
  protected:

    Cluster_Definitions_Base * p_cluster;
    
    std::string m_name;
    ATOOLS::Event_Weights m_weights;
    int m_on, m_kttype;

  public:

    Shower_Base(const std::string &name);

    virtual ~Shower_Base();

    virtual int  PerformShowers(ATOOLS::Cluster_Amplitude *const ampl) = 0;
    virtual void CleanUp()=0;

    virtual Cluster_Definitions_Base *GetClusterDefinitions() = 0;

    virtual bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl,
                               const bool & soft=false) = 0;

    virtual void SetRemnants(REMNANTS::Remnant_Handler * remnants) {}

    static void ShowSyntax(const int mode);

    inline const std::string &Name() const { return m_name; }

    inline ATOOLS::Event_Weights Weights() const { return m_weights; }

    inline int On() { return m_on; }

    inline void SetOn(const int on) { m_on=on; }

    inline int KTType() const { return m_kttype; }

  };// end of class Shower_Base

  struct Shower_Key {
    MODEL::Model_Base *p_model;
    ISR_Handler *p_isr;
    inline Shower_Key(MODEL::Model_Base *const model,
		      ISR_Handler *const isr):
      p_model(model), p_isr(isr) {}
  };//end of struct Shower_Key

  typedef ATOOLS::Getter_Function
  <Shower_Base,Shower_Key> Shower_Getter;

}// end of namespace ATOOLS

#endif
