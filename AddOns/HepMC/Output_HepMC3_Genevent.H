#ifndef SHERPA_Tools_Output_HepMC3_Genevent_H
#define SHERPA_Tools_Output_HepMC3_Genevent_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepMC3_Interface.H"
#include <iostream>
#include <memory>

namespace HepMC3 {
  class GenEvent;
  class GenVertex;
  class GenParticle;
  class GenCrossSection;
  class Writer;
  class WriterAscii;
  class WriterAsciiHepMC2;
  class WriterHEPEVT;
#ifdef USING__HEPMC3__WRITERROOTTREE
  class WriterRootTree;
#endif
#ifdef USING__HEPMC3__WRITERROOT
  class WriterRoot;
#endif
}
namespace HepMC = HepMC3;

namespace SHERPA {

  class Output_HepMC3_Genevent : public Output_Base {
    HepMC3_Interface m_hepmc3;
    HepMC::Writer *p_writer=0;
    std::shared_ptr<HepMC::GenCrossSection> p_xs;
    std::shared_ptr<HepMC::GenRunInfo> m_run_info;
    std::string m_basename;
    int m_iotype;
  public:
    Output_HepMC3_Genevent(const Output_Arguments &args);
    ~Output_HepMC3_Genevent();
    void SetXS(const double& xs, const double& xserr);
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile();
  };
}
#endif
