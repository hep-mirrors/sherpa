#ifndef SHERPA_Tools_Output_HepEvt_H
#define SHERPA_Tools_Output_HepEvt_H

#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepEvt_Interface.H"

namespace SHERPA {

  class Output_HepEvt : public Output_Base {
    HepEvt_Interface m_hepevt;
  public:
    Output_HepEvt(std::string basename,std::string ext,int precision) :
      Output_Base(basename, ext, precision) {}
    void Output(ATOOLS::Blob_List* blobs, double weight);
  };

}

#endif
