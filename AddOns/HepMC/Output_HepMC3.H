#ifndef SHERPA_Tools_Output_HepMC3_H
#define SHERPA_Tools_Output_HepMC3_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepMC3_Interface.H"
#include <iostream>
#include <memory>
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#endif

namespace HepMC3 {
  class GenRunInfo;
  class GenEvent;
  class GenVertex;
  class GenParticle;
  class Writer;
}
namespace HepMC = HepMC3;

namespace SHERPA {

  class Output_HepMC3 : public Output_Base {
    HepMC3_Interface m_hepmc3;
    HepMC::Writer* p_writer = nullptr;
    std::shared_ptr<HepMC::GenRunInfo> m_run_info;
    std::string m_basename, m_ext;
#ifdef USING__GZIP
    ATOOLS::ogzstream m_outstream;
#else
    std::ofstream     m_outstream;
#endif
    int m_iotype;
    ATOOLS::Weights_Map m_xs;
    ATOOLS::Weights_Map m_err;

  protected:
    bool m_short;

  public:
    Output_HepMC3(const Output_Arguments &args);
    ~Output_HepMC3();
    void SetXS(const ATOOLS::Weights_Map& xs,
               const ATOOLS::Weights_Map& err) override;
    void Output(ATOOLS::Blob_List* blobs) override;
    void ChangeFile() override;
  };

  class Output_HepMC3_GenEvent : public Output_HepMC3 {
  public:
    Output_HepMC3_GenEvent(const Output_Arguments &args) : Output_HepMC3(args) {
      m_short = false;
    }
  };

  class Output_HepMC3_Short : public Output_HepMC3 {
  public:
    Output_HepMC3_Short(const Output_Arguments &args) : Output_HepMC3(args) {
      m_short = true;
    }
  };

}

#endif
