#ifndef PHASIC_Process_ME_Generator_Base_H
#define PHASIC_Process_ME_Generator_Base_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS { class Cluster_Amplitude; }
namespace BEAM { class Beam_Spectra_Handler; }
namespace PDF { class ISR_Handler; }
namespace YFS { class YFS_Handler; }
namespace MODEL { class Model_Base; }
namespace REMNANTS { class Remnant_Handler; }

namespace PHASIC {

  class ME_Generators;

  class ME_Generator_Base: public ATOOLS::Mass_Selector
  {
  private:

    std::string         m_name;
    ATOOLS::Flavour_Set m_psmass;
    ME_Generators*      p_gens;
    REMNANTS::Remnant_Handler *p_remnant;

    void RegisterDefaults();
    void RegisterDipoleParameters();
    void RegisterNLOParameters();
    template <typename T>
    void SetParameter(const std::string& param, const T& def);

  protected:

    int  m_massmode;
    void SetPSMasses();

  public:

    // constructor
    ME_Generator_Base(const std::string &name);

    // destructor
    virtual ~ME_Generator_Base();

    // member functions
    virtual bool Initialize(MODEL::Model_Base* const,
			    BEAM::Beam_Spectra_Handler* const,
			    PDF::ISR_Handler* const,
			    YFS::YFS_Handler* const) = 0;

    virtual Process_Base *InitializeProcess(const Process_Info &pi,
                                            bool add)=0;
    virtual int PerformTests() = 0;
    virtual bool NewLibraries() = 0;

    Process_Base *InitializeProcess
    (ATOOLS::Cluster_Amplitude *const ampl,const int mode,
     const std::string &gen="",const std::string &addname="");

    int ShiftMasses(ATOOLS::Cluster_Amplitude *const ampl);

    double Mass(const ATOOLS::Flavour &fl) const;

    static void ShowSyntax(const int mode);

    // inline functions
    inline std::string Name() const { return m_name; }

    inline int SetMassMode(const int mode)
    { int mm(m_massmode); m_massmode=mode; return mm; }

    inline int MassMode() { return m_massmode; }

    inline void SetGenerators(ME_Generators *const gens) { p_gens=gens; }

    inline ME_Generators *Generators() const { return p_gens; }

    virtual void SetModel(MODEL::Model_Base*) {};

    inline void SetRemnantHandler(REMNANTS::Remnant_Handler *remnant) {
      p_remnant = remnant;
    }
  };// end of class ME_Generator_Base

  struct ME_Generator_Key {
    ME_Generator_Key() {};
  };//end of struct ME_Generator_Key

  typedef ATOOLS::Getter_Function
  <ME_Generator_Base,ME_Generator_Key> ME_Generator_Getter;

}// end of namespace PHASIC

#endif
