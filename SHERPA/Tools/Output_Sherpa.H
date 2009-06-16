#ifndef SHERPA_Tools_Output_Sherpa_H
#define SHERPA_Tools_Output_Sherpa_H

#include "SHERPA/Tools/Output_Base.H"

namespace SHERPA {

  class Output_Sherpa : public Output_Base {
  public:
    Output_Sherpa(std::string basename,std::string ext,int precision) :
      Output_Base(basename, ext, precision) {}
    void Header();
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void Footer(std::string nr);
  };

}

#endif
