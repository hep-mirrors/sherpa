#ifndef SHERPA_Tools_Output_Base_H
#define SHERPA_Tools_Output_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Weights.H"

#include <string>

namespace SHERPA {

  struct Output_Arguments {
    std::string m_outpath, m_outfile;
    Output_Arguments(const std::string outpath,
		     const std::string outfile):
      m_outpath(outpath), m_outfile(outfile)
    {}
  };

  class Event_Handler;

  class Output_Base {
  protected:

    std::string m_name;

    Event_Handler *p_eventhandler;

  public:

    typedef ATOOLS::Getter_Function
    <Output_Base,Output_Arguments> Getter_Function;

  public:

    Output_Base(const std::string &name);

    virtual ~Output_Base();

    virtual void Header();
    virtual void Footer();

    virtual void ChangeFile();

    virtual void SetXS(const ATOOLS::Weights_Map& xs,
		       const ATOOLS::Weights_Map& err);
    virtual void Output(ATOOLS::Blob_List* blobs) = 0;

    inline const std::string Name() { return m_name; }

    void SetEventHandler(Event_Handler *const eh) { p_eventhandler = eh; }

  };

  typedef std::vector<Output_Base*> Output_Vector;

}

#endif
