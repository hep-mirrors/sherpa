#ifndef SHERPA_Tools_Event_Reader_H
#define SHERPA_Tools_Event_Reader_H

#include <fstream>
#include "SHERPA/Tools/Event_Reader_Base.H"


namespace SHERPA {
  class Event_Reader : public Event_Reader_Base {
    std::ifstream*     p_instream;

    std::string        m_generator;
    bool               f_gz;
    int                m_inputmode, m_eventmode, m_phasemode;
    int                m_evtnumber;

    void               InitialSettings();
    void               CloseFile();
    bool               ReadInEvent(ATOOLS::Blob_List *);
    bool               ReadInSimpleHepEvtEvent(ATOOLS::Blob_List *);
    ATOOLS::Particle * TranslateFromInput(std::string,int&,const int=0);
  public:
    Event_Reader(const std::string & =std::string(""),const std::string & =std::string(""));
    bool        FillBlobs(ATOOLS::Blob_List *);
  };
}

#endif
