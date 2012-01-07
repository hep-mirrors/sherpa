#ifndef SHERPA_Tools_RootNtuple_Reader_H
#define SHERPA_Tools_RootNtuple_Reader_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Event_Reader_Base.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include <vector>

namespace PDF { class ISR_Handler; }

namespace MODEL { class Model_Base; }

namespace PHASIC { class Scale_Setter_Base; }

namespace SHERPA {

  struct RootNTupleReader_Variables;

  class RootNtuple_Reader : public Event_Reader_Base {
    RootNTupleReader_Variables *p_vars;
    int                m_eventmode,m_check,m_calc,m_oqcd;
    size_t             m_evtid,m_subevtid,m_evtcnt,m_entries,m_evtpos;
    ATOOLS::NLO_subevtlist m_nlos;
    PHASIC::Scale_Setter_Arguments m_sargs;
    std::map<int,PHASIC::Process_Base*> m_procs;
    std::map<int,PHASIC::Scale_Setter_Base*> m_scales;
    PDF::ISR_Handler *p_isr;

    void               CloseFile();
    bool               ReadInEntry();
    bool               ReadInFullEvent(ATOOLS::Blob_List *);
    bool               ReadInSubEvent(ATOOLS::Blob_List *);
    double CalculateWeight(const double &mur2,const double &muf2,const int mode) const;
  public:
    RootNtuple_Reader(const std::string &,const std::string &,int mode,
		      MODEL::Model_Base *const model,PDF::ISR_Handler *const isr);
    ~RootNtuple_Reader();
    bool        FillBlobs(ATOOLS::Blob_List *);
  };
}

#endif
