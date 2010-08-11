#ifndef SHERPA_Tools_Analysis_Interface_H
#define SHERPA_Tools_Analysis_Interface_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace SHERPA {
  
  class Event_Handler;

  struct Analysis_Arguments {
    std::string m_inpath, m_infile, m_outpath;
    Analysis_Arguments(const std::string inpath,
		       const std::string infile,
		       const std::string outpath):
      m_inpath(inpath), m_infile(infile), 
      m_outpath(outpath) {}
  };// end of struct Analysis_Arguments

  class Analysis_Interface {
  protected:
    std::string m_name;
    const Event_Handler* p_eventhandler;
  public:

    typedef ATOOLS::Getter_Function
    <Analysis_Interface,Analysis_Arguments> 
    Analysis_Getter_Function;

  public:

    inline Analysis_Interface(const std::string &name):
      m_name(name) {}

    virtual ~Analysis_Interface();

    virtual bool Init() = 0;
    virtual bool Run(ATOOLS::Blob_List *const bl) = 0;
    virtual bool Finish() = 0;

    virtual void CleanUp();
    virtual bool WriteOut();

    virtual void ShowSyntax(const int mode);
    
    void SetEventHandler(const Event_Handler* eh) { p_eventhandler = eh; }

    inline std::string Name() const { return m_name; }

  };// end of class Analysis_Interface

  typedef std::vector<Analysis_Interface*> Analysis_Vector;

}// end of namespace SHERPA

#endif
