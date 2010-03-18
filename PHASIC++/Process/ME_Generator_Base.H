#ifndef PHASIC_Process_ME_Generator_Base_H
#define PHASIC_Process_ME_Generator_Base_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS { class Cluster_Amplitude; }
namespace BEAM { class Beam_Spectra_Handler; }
namespace PDF {
  class Cluster_Definitions_Base; 
  class ISR_Handler;
}
namespace MODEL { class Model_Base; }

namespace PHASIC {

  class ME_Generator_Base: public ATOOLS::Mass_Selector {
  private:

    std::string m_name;

  public :

    // constructor
    inline ME_Generator_Base(const std::string &name): m_name(name) {}

    // destructor
    virtual ~ME_Generator_Base();

    // member functions
    virtual bool Initialize(const std::string &path,const std::string &file,
			    MODEL::Model_Base *const model,
			    BEAM::Beam_Spectra_Handler *const beam,
			    PDF::ISR_Handler *const isr) = 0;
    virtual Process_Base *InitializeProcess(const Process_Info &pi, bool add)=0;
    virtual bool PerformTests() = 0;

    virtual Process_Base *GetProcesses() = 0;

    virtual void SetClusterDefinitions
    (PDF::Cluster_Definitions_Base *const defs) = 0;

    virtual ATOOLS::Cluster_Amplitude *
    ClusterConfiguration(Process_Base *const proc,
			 const size_t &mode=0,const double &kt2=0.0) = 0;

    double Mass(const ATOOLS::Flavour &fl) const;

    static void ShowSyntax(const int mode);

    // inline functions
    inline std::string Name() const { return m_name; }

  };// end of class ME_Generator_Base

  struct ME_Generator_Key {
  };//end of struct ME_Generator_Key

  typedef ATOOLS::Getter_Function
  <ME_Generator_Base,ME_Generator_Key> ME_Generator_Getter;

}// end of namespace PHASIC

#endif
