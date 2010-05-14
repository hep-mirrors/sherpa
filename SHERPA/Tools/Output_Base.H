#ifndef SHERPA_Tools_Output_Base_H
#define SHERPA_Tools_Output_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include <string>
#include <fstream>
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#endif

namespace SHERPA {

  class Output_Base {
  protected:
    std::string m_basename, m_ext;
#ifdef USING__GZIP
    ATOOLS::ogzstream m_outstream;
#else
    std::ofstream     m_outstream;
#endif

  public:
    Output_Base() {}
    Output_Base(std::string basename,std::string ext,int precision);
    virtual ~Output_Base();

    virtual void ChangeFile(std::string nr);

    inline const std::string& BaseName() const { return m_basename; }
    inline const std::string& Ext() const { return m_ext; }

    virtual void Header() {}
    virtual void Footer(std::string) {}
    virtual void SetXS(const double&, const double&) {}
    virtual void Output(ATOOLS::Blob_List* blobs, double weight) = 0;
  };
}

#endif
