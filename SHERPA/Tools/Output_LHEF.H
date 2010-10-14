#ifndef Output_LHEF_H
#define Output_LHEF_H

#include "Output_Base.H"

namespace SHERPA {

  class Output_LHEF : public Output_Base {
    double m_xs,m_xserr,m_max;
    
  public:
    Output_LHEF(std::string basename,std::string ext,int precision) :
      Output_Base(basename, ext, precision) {}
    void SetXS(const double& xs, const double& xserr);
    void Header();
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void Footer(std::string nr);
  };

}

#endif
