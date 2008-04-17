#ifndef Output_HepMC2_Ascii_H
#define Output_HepMC2_Ascii_H

#include "Output_Base.H"
#include "HepMC2_Interface.H"
#include "HepMC/IO_Ascii.h"

namespace SHERPA {

  class Output_HepMC2_Ascii : public Output_Base {
    HepMC2_Interface m_hepmc2;
    HepMC::IO_Ascii* p_ioascii;
  public:
    Output_HepMC2_Ascii(std::string basename,std::string ext,int precision);
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile(std::string number);
  };

}

#endif
