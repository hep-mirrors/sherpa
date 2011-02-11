#ifndef SHERPA_Tools_RootNtuple_Reader_H
#define SHERPA_Tools_RootNtuple_Reader_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Event_Reader_Base.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include <vector>

namespace SHERPA {

  struct RootNTupleReader_Variables;

  class RootNtuple_Reader : public Event_Reader_Base {
    RootNTupleReader_Variables *p_vars;
    int                m_eventmode;
    size_t             m_evtid,m_subevtid,m_evtcnt,m_entries,m_evtpos;
    ATOOLS::NLO_subevtlist m_nlos;

    void               CloseFile();
    bool               ReadInEntry();
    bool               ReadInFullEvent(ATOOLS::Blob_List *);
    bool               ReadInSubEvent(ATOOLS::Blob_List *);
  public:
    RootNtuple_Reader(const std::string & =std::string(""),const std::string & =std::string(""),int mode=0);
    bool        FillBlobs(ATOOLS::Blob_List *);
  };
}

#endif
