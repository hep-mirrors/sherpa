#ifndef SHERPA_Tools_Output_HepEvt_H
#define SHERPA_Tools_Output_HepEvt_H

#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepEvt_Interface.H"
#include <fstream>
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#endif

namespace SHERPA {

  class Output_HepEvt : public Output_Base {
    HepEvt_Interface m_hepevt;
    std::string m_basename, m_ext;
#ifdef USING__GZIP
    ATOOLS::ogzstream m_outstream;
#else
    std::ofstream     m_outstream;
#endif    
  public:
    Output_HepEvt(const Output_Arguments &args);
    ~Output_HepEvt();
    void ChangeFile();
    void Output(ATOOLS::Blob_List* blobs, double weight);
  };

}

#endif
