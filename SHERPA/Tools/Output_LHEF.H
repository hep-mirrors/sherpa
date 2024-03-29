#ifndef SHERPA_Tools_Output_LHEF_H
#define SHERPA_Tools_Output_LHEF_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#include <fstream>
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#endif

namespace SHERPA {

  class Output_LHEF : public Output_Base {
    double m_xs,m_xserr,m_max;
    std::string m_basename, m_ext;
#ifdef USING__GZIP
    ATOOLS::ogzstream m_outstream;
#else
    std::ofstream     m_outstream;
#endif    
    int m_bntp;
  public:
    Output_LHEF(const Output_Arguments &args);
    ~Output_LHEF();
    void Header() override;
    void Footer() override;
    void ChangeFile() override;
    void SetXS(const ATOOLS::Weights_Map& xs,
	       const ATOOLS::Weights_Map& err) override;
    void Output(ATOOLS::Blob_List* blobs) override;
  private:
    void RegisterDefaults() const;
  };

}

#endif
