#ifndef SHERPA_Tools_Output_HepMC2_Short_H
#define SHERPA_Tools_Output_HepMC2_Short_H

#include "SHERPA/Tools/Output_Base.H"
#include "AddOns/HepMC/HepMC2_Interface.H"
#include <fstream>
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#endif

namespace HepMC { class IO_GenEvent; class GenCrossSection; class GenEvent; }

namespace SHERPA {

  class Output_HepMC2_Short : public Output_Base {
    HepMC2_Interface m_hepmc2;
    HepMC::IO_GenEvent* p_iogenevent;
    HepMC::GenCrossSection* p_xs;
    HepMC::GenEvent* p_event;
    std::string m_basename, m_ext;
#ifdef USING__GZIP
    ATOOLS::ogzstream m_outstream;
#else
    std::ofstream     m_outstream;
#endif
  public:
    Output_HepMC2_Short(const Output_Arguments &args);
    ~Output_HepMC2_Short();
    void SetXS(const double& xs, const double& xserr);
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile();
  };

}

#endif
