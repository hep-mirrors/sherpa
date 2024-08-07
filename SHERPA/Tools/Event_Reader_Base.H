#ifndef SHERPA_Tools_Event_Reader_Base_H
#define SHERPA_Tools_Event_Reader_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS { class Variations; }

namespace MODEL { class Model_Base; }

namespace PDF { class ISR_Handler; }

namespace YFS { class YFS_Handler; }

namespace SHERPA {

  struct Input_Arguments {
    std::string m_inpath, m_infile;
    MODEL::Model_Base *p_model;
    PDF::ISR_Handler *p_isr;
    YFS::YFS_Handler *p_yfs;
    Input_Arguments(const std::string inpath,
		    const std::string infile,
		    MODEL::Model_Base *const model,
		    PDF::ISR_Handler *const isr,
        YFS::YFS_Handler *const yfs):
      m_inpath(inpath), m_infile(infile),
      p_model(model), p_isr(isr), p_yfs(yfs) {}
  };// end of struct Input_Arguments

  class Event_Reader_Base {
  public:
    typedef ATOOLS::Getter_Function
    <Event_Reader_Base,Input_Arguments> Getter_Function;
  protected:
    std::string        m_path, m_file;

    ATOOLS::Variations * p_variations;

  public:
    Event_Reader_Base(const Input_Arguments &args)
      : m_path(args.m_inpath), m_file(args.m_infile), p_variations(NULL) {}
    virtual ~Event_Reader_Base();
    virtual bool FillBlobs(ATOOLS::Blob_List *)=0;

    std::string GetPath() { return m_path; }

    void        SetVariations(ATOOLS::Variations * v) { p_variations = v; }
  };
}

#endif
