#ifndef SHERPA_Tools_Output_LHEF_H
#define SHERPA_Tools_Output_LHEF_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#include <fstream>
#include "ATOOLS/Org/Gzip_Stream.H"

namespace SHERPA {

  class Output_LHEF : public Output_Base {
    double m_xs,m_xserr,m_max;
    std::string m_basename, m_ext;
    ATOOLS::Gzip_Stream m_outstream;
    int m_bntp;
  public:
    Output_LHEF(const Output_Arguments &args);
    ~Output_LHEF();
    void ChangeFile();
    void SetXS(const double& xs, const double& xserr);
    void Header();
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void Footer();
  };

}

#endif
