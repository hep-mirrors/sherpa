#ifndef SHERPA_Tools_Event_Reader_Base_H
#define SHERPA_Tools_Event_Reader_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS { class Data_Reader; }

namespace MODEL { class Model_Base; }

namespace PDF { class ISR_Handler; }

namespace SHERPA {

  struct Input_Arguments {
    std::string m_inpath, m_infile;
    MODEL::Model_Base *p_model;
    PDF::ISR_Handler *p_isr;
    ATOOLS::Data_Reader *p_reader;
    Input_Arguments(const std::string inpath,
		    const std::string infile,
		    ATOOLS::Data_Reader *const read,
		    MODEL::Model_Base *const model,
		    PDF::ISR_Handler *const isr):
      m_inpath(inpath), m_infile(infile),
      p_reader(read), p_model(model), p_isr(isr) {}
  };// end of struct Input_Arguments

  class Event_Reader_Base {
  public:
    typedef ATOOLS::Getter_Function
    <Event_Reader_Base,Input_Arguments> Getter_Function;
  protected:
    std::string        m_path, m_file;

    double             m_weight;

  public:
    Event_Reader_Base(const Input_Arguments &args)
      : m_path(args.m_inpath), m_file(args.m_infile) {}
    virtual ~Event_Reader_Base();
    virtual bool FillBlobs(ATOOLS::Blob_List *)=0;

    double      Weight()  { return m_weight; }
    std::string GetPath() { return m_path; }
  };
}

#endif
