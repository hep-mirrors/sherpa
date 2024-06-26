#ifndef PHASIC_Process_ME_Generators_H
#define PHASIC_Process_ME_Generators_H

#include <vector>
#include <string>

namespace MODEL    { class Model_Base; }
namespace BEAM     { class Beam_Spectra_Handler; }
namespace PDF      { class ISR_Handler; }
namespace YFS      { class YFS_Handler; }
namespace REMNANTS { class Remnant_Handler; }

namespace PHASIC   {

  class Process_Base;
  class Process_Info;
  class ME_Generator_Base;

  class ME_Generators : public std::vector<ME_Generator_Base*> {
  private:

    MODEL::Model_Base *p_model;

    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;
    YFS::YFS_Handler           *p_yfs;

  public:

    ME_Generators();
    ~ME_Generators();

    bool InitializeGenerators(MODEL::Model_Base *model,
                              BEAM::Beam_Spectra_Handler *beam,
                              PDF::ISR_Handler *isr,
                              YFS::YFS_Handler *yfs);
    int PerformTests();
    bool NewLibraries();

    bool LoadGenerator(const std::string &name);

    Process_Base* InitializeProcess(const Process_Info &pi, bool add);

    inline MODEL::Model_Base *Model() const { return p_model; }
    void SetModel(MODEL::Model_Base*);
    void SetRemnant(REMNANTS::Remnant_Handler *remnant);
  };

}

#endif
