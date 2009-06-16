#ifndef PDF_Main_Shower_Base_H
#define PDF_Main_Shower_Base_H

#include "ATOOLS/Org/Getter_Function.H"

#include <string>

namespace ATOOLS {

  class Cluster_Definitions_Base;
  class Cluster_Amplitude;
  class Blob_List;
  class Data_Reader;

}

namespace MODEL { class Model_Base; }

namespace PDF {

  class ISR_Handler;

  class Shower_Base {
  protected:

    ATOOLS::Cluster_Definitions_Base * p_cluster;
    
    std::string m_name;

  public:

    Shower_Base(const std::string &name);

    virtual ~Shower_Base();

    virtual int  PerformShowers() = 0;
    virtual int  PerformDecayShowers() = 0;
    virtual int  TrialEmission() = 0;

    virtual double CouplingWeight(ATOOLS::Cluster_Amplitude *const ampl) = 0;
    virtual double TrialWeight(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    virtual bool ExtractPartons(ATOOLS::Blob_List *const bl) = 0;
    virtual void CleanUp()=0;

    virtual ATOOLS::Cluster_Definitions_Base *GetClusterDefinitions() = 0;

    virtual ATOOLS::Cluster_Amplitude *GetRealEmissionAmplitude() = 0;

    virtual bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    virtual double CalculateWeight(ATOOLS::Cluster_Amplitude *const ampl) = 0;

    virtual std::string GetKT2(const std::string &jm2) const = 0;

    static void ShowSyntax(const int mode);

    inline const std::string &Name() const { return m_name; }

  };// end of class Shower_Base

  struct Shower_Key {
    MODEL::Model_Base *p_model;
    ISR_Handler *p_isr;
    ATOOLS::Data_Reader *p_read;
    inline Shower_Key(MODEL::Model_Base *const model,
		      ISR_Handler *const isr,
		      ATOOLS::Data_Reader *const read):
      p_model(model), p_isr(isr), p_read(read) {}
  };//end of struct Shower_Key

  typedef ATOOLS::Getter_Function
  <Shower_Base,Shower_Key> Shower_Getter;

}// end of namespace ATOOLS

#endif
