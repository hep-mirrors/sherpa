#ifndef EvtReadin_Phase_H
#define EvtReadin_Phase_H

#include "Blob_List.H"
#include "Event_Reader.H"
#include "Output_Handler.H"
#include "Event_Phase_Handler.H"

namespace SHERPA {
  class EvtReadin_Phase : public Event_Phase_Handler {
    IO_Handler     * p_iohandler;
    Event_Reader   * p_evtreader;
    std::string      m_path;
    bool             m_read;
  public:
    EvtReadin_Phase(Event_Reader *);
    EvtReadin_Phase(IO_Handler *);
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string &);
  };
}

#endif
