#ifndef SHERPA_Tools_Output_Base_H
#define SHERPA_Tools_Output_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Getter_Function.H"
#include <string>

namespace ATOOLS { class Data_Reader; }

namespace SHERPA {

  struct Output_Arguments {
    std::string m_outpath, m_outfile;
    ATOOLS::Data_Reader *p_reader;
    Output_Arguments(const std::string outpath,
		     const std::string outfile,
		     ATOOLS::Data_Reader *const read):
      m_outpath(outpath), m_outfile(outfile), p_reader(read) {}
  };// end of struct Output_Arguments

  class Event_Handler;

  class Output_Base {
  protected:

    std::string m_name;

    const Event_Handler* p_eventhandler;

  public:

    typedef ATOOLS::Getter_Function
    <Output_Base,Output_Arguments> Getter_Function;

  public:

    Output_Base(const std::string &name);

    virtual ~Output_Base();

    virtual void Header();
    virtual void Footer();

    virtual void ChangeFile();

    virtual void SetXS(const double&, const double&);
    virtual void Output(ATOOLS::Blob_List* blobs, double weight) = 0;

    inline const std::string Name() { return m_name; }

    void SetEventHandler(const Event_Handler* eh) { p_eventhandler = eh; }

  };

  typedef std::vector<Output_Base*> Output_Vector;

}

#endif
