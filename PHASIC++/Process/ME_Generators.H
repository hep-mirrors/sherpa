#ifndef PHASIC_Process_ME_Generators_H
#define PHASIC_Process_ME_Generators_H

#include <vector>
#include <string>

namespace MODEL    { class Model_Base; }
namespace BEAM     { class Beam_Spectra_Handler; }
namespace PDF      { class ISR_Handler; }

namespace PHASIC   {
  
  class Process_Base;
  class Process_Info;
  class ME_Generator_Base;

  class ME_Generators : public std::vector<ME_Generator_Base*> {
  private:

    std::string m_path, m_file;
    MODEL::Model_Base *p_model;

  public:

    ME_Generators(const std::string &path,const std::string &file);
    ~ME_Generators();

    bool InitializeGenerators(MODEL::Model_Base *model,
                              BEAM::Beam_Spectra_Handler *beam,
                              PDF::ISR_Handler *isr);
    int PerformTests();
    bool NewLibraries();
    Process_Base* InitializeProcess(const Process_Info &pi, bool add);

    inline MODEL::Model_Base *Model() const { return p_model; }

  };

}

#endif
