#ifndef PHASIC_Process_ME_Generator_Base_H
#define PHASIC_Process_ME_Generator_Base_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS { class Cluster_Amplitude; }
namespace BEAM { class Beam_Spectra_Handler; }
namespace PDF { class ISR_Handler; }
namespace YFS { class YFS_Handler; }
namespace MODEL { class Model_Base; }

namespace PHASIC {

  class ME_Generators;

  /// Momentum mapping that puts a 2->n amplitude on the parton-shower mass
  /// shell, with the incoming legs it treats specially:
  /// - dis: m_fixed is the lepton of a lepton beam, m_pdf the parton it hits;
  /// - epa: m_fixed is a non-lepton without PDF (e.g. a real photon from an
  ///   EPA spectrum), m_pdf the leg from a PDF it hits;
  /// - standard: both legs are rescaled alike, m_fixed = m_pdf = -1.
  struct Mass_Shift_Mapping {
    enum code { standard, dis, epa };
    code m_mode;
    int m_fixed, m_pdf;
  };

  /// Select the mapping from the incoming flavours, the beams with a PDF and
  /// the beams whose bunch is a lepton (bit i set for beam i in both masks).
  Mass_Shift_Mapping SelectMassShiftMapping(const ATOOLS::Flavour &in0,
                                            const ATOOLS::Flavour &in1,
                                            const int pdfs,
                                            const int leptonbunches);

  class ME_Generator_Base: public ATOOLS::Mass_Selector
  {
  private:

    std::string         m_name;
    ATOOLS::Flavour_Set m_psmass;
    ME_Generators*      p_gens;

    void RegisterDefaults();
    void RegisterDipoleParameters();
    void RegisterNLOParameters();
    template <typename T>
    void SetParameter(const std::string& param, const T& def);

    int ShiftMassesDefault(ATOOLS::Cluster_Amplitude *const ampl,
                           ATOOLS::Vec4D cms, const int pdfs);
    int ShiftMassesDIS(ATOOLS::Cluster_Amplitude *const ampl,
                       const int lep_leg, const int had_leg);
    int ShiftMassesEPA(ATOOLS::Cluster_Amplitude *const ampl,
                       const int lep_leg, const int had_leg);

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
  };// end of class ME_Generator_Base

  struct ME_Generator_Key {
    ME_Generator_Key() {};
  };//end of struct ME_Generator_Key

  typedef ATOOLS::Getter_Function
  <ME_Generator_Base,ME_Generator_Key> ME_Generator_Getter;

}// end of namespace PHASIC

#endif
