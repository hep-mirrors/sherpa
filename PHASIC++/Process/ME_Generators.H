#ifndef ME_GENERATORS_H
#define ME_GENERATORS_H

#include <vector>
#include <string>

namespace MODEL    { class Model_Base; }
namespace BEAM     { class Beam_Spectra_Handler; }
namespace PDF      { class ISR_Handler; }

namespace PHASIC   {
  
  class Process_Base;
  class Process_Info;
  class ME_Generator_Base;

  class ME_Generators : public std::vector<ME_Generator_Base*>
  {
    std::string m_path, m_file;
  public:
    ME_Generators(const std::string &path,const std::string &file);
    ~ME_Generators();

    bool InitializeGenerators(MODEL::Model_Base *model,
                              BEAM::Beam_Spectra_Handler *beam,
                              PDF::ISR_Handler *isr);
    bool PerformTests();
    Process_Base* InitializeProcess(const Process_Info &pi, bool add);
  };

}

#endif // ME_GENERATORS_H
