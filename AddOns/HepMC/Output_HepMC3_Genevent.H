#ifndef SHERPA_Tools_Output_HepMC3_Genevent_H
#define SHERPA_Tools_Output_HepMC3_Genevent_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepMC3_Interface.H"

namespace HepMC {
  class GenEvent;
  class GenVertex;
  class GenParticle;
  class GenCrossSection;
  class Writer;
  class WriterAscii;
  class WriterHEPEVT;
#ifdef USING__HEPMC3__WRITERROOTTREE
  class WriterRootTree;
#endif
#ifdef USING__HEPMC3__WRITERROOT
  class WriterRoot;
#endif
}


namespace SHERPA {

  class Output_HepMC3_Genevent : public Output_Base {
    HepMC3_Interface m_hepmc3;
    HepMC::Writer *p_writer=0;
    std::string m_basename, m_ext;
    int m_iotype;
  public:
    Output_HepMC3_Genevent(const Output_Arguments &args);
    ~Output_HepMC3_Genevent();
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile();
  };
}
#endif
