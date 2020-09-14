#ifndef SHERPA_Single_Events_Userhook_Phase_H
#define SHERPA_Single_Events_Userhook_Phase_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/Tools/Userhook_Base.H"

namespace SHERPA {

  class Sherpa;

  class Userhook_Phase: public Event_Phase_Handler {
  private:
    
    Userhook_Vector m_userhooks;

  public:
    
    // constructor
    Userhook_Phase(Sherpa* sherpa);
    ~Userhook_Phase();

    // member functions
    void InitializeHooks(Sherpa* sherpa);
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void Finish(const std::string &);
    void CleanUp(const size_t & mode=0);

  };// end of class Userhook_Phase

}// end of namespace SHERPA

#endif
