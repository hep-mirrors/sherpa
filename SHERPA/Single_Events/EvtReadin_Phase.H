#ifndef SHERPA_Single_Events_EvtReadin_Phase_H
#define SHERPA_Single_Events_EvtReadin_Phase_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Variations.H"
#include "SHERPA/Tools/Event_Reader_Base.H"
#include "SHERPA/Single_Events/Event_Phase_Handler.H"

namespace SHERPA {

  class EvtReadin_Phase : public Event_Phase_Handler {
    Event_Reader_Base    * p_evtreader;
    std::string            m_path;
    bool                   m_read;
  public:
    EvtReadin_Phase(Event_Reader_Base *, ATOOLS::Variations *);
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };
}

#endif
