#ifndef SHERPA_Tools_Analysis_Interface_H
#define SHERPA_Tools_Analysis_Interface_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace SHERPA {

  struct Analysis_Arguments {
    std::string m_inpath, m_infile, m_outpath;
    Analysis_Arguments(const std::string inpath,
		       const std::string infile,
		       const std::string outpath):
      m_inpath(inpath), m_infile(infile), 
      m_outpath(outpath) {}
  };// end of struct Analysis_Arguments

  class Analysis_Interface {
  public:

    typedef ATOOLS::Getter_Function
    <Analysis_Interface,Analysis_Arguments> 
    Analysis_Getter_Function;

  public:

    virtual ~Analysis_Interface();

    virtual bool Init() = 0;
    virtual bool Run(ATOOLS::Blob_List *const bl) = 0;
    virtual bool Finish() = 0;

    virtual void CleanUp();
    virtual bool WriteOut();

    virtual void ShowSyntax(const int mode);

  };// end of class Analysis_Interface

}// end of namespace SHERPA

#endif
