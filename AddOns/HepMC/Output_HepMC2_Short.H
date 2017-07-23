#ifndef SHERPA_Tools_Output_HepMC2_Short_H
#define SHERPA_Tools_Output_HepMC2_Short_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepMC2_Interface.H"
#include <fstream>
#include "ATOOLS/Org/Gzip_Stream.H"

namespace HepMC { class IO_GenEvent; class GenCrossSection; class GenEvent; }

namespace SHERPA {

  class Output_HepMC2_Short : public Output_Base {
    HepMC2_Interface m_hepmc2;
    HepMC::IO_GenEvent* p_iogenevent;
    HepMC::GenCrossSection* p_xs;
    HepMC::GenEvent* p_event;
    std::string m_basename, m_ext;
    ATOOLS::Gzip_Stream m_outstream;
  public:
    Output_HepMC2_Short(const Output_Arguments &args);
    ~Output_HepMC2_Short();
    void SetXS(const double& xs, const double& xserr);
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile();
  };

}

#endif
