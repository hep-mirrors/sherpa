#ifndef Event_Phase_Handler_h
#define Event_Phase_Handler_h

//#include "Parton_List.H"
//#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Event_Phase_Handler {
  protected :
    std::string m_type,m_name;
  public :
    Event_Phase_Handler();
    virtual ~Event_Phase_Handler();
    virtual bool Treat(ATOOLS::Blob_List *, double &) = 0;
    virtual void        CleanUp()                        = 0;
    virtual std::string Type();
    virtual void        SetType(std::string);
    virtual std::string Name();
    virtual void        SetName(std::string);
  };

  /*!
    \class Event_Phase_Handler
    \brief The base class for all specific event phases (such as, e.g., JetEvolution)

    This is the - mainly abstract - base class for all event phases. Apart from the name and type and 
    methods to access these characteristics, and apart from constructors and destructors, this class 
    contains two central methods, namely Treat(ATOOLS::Blob_List *, double &) and CleanUp(). The 
    former is the one through which the concrete event phases interact during the event generation. This 
    happens by modifying individual blobs of the list or adding or removing blobs from the list. The
    latter method is used to reset the event phase handlers between two events.
  */
  /*!
    \var std::string Event_Phase_Handler::m_type
    The type of the phase. At the moment two types are implemented, "Perturbative" and "Hadronization".
  */
  /*!
    \var std::string Event_Phase_Handler::m_name
    The name of the phase. The combination of type and name helps to ensure that only one Event_Phase_Handler
    of each type enters into the Event_Handler.
  */
  /*!
    \fn Event_Phase_Handler::Event_Phase_Handler()
    Here, the type and name are set - by default - to "Unspecified" and "No Name", respectively.
  */
  /*!
    \fn virtual bool Event_Phase_Handler::Treat(ATOOLS::Blob_List *, double &)
    Virtual method to be filled by the specific Event_Phase_Handler. this method basically defines
    the role, any of these handlers plays in the event generation process. 
  */
  /*!
    \fn virtual void Event_Phase_Handler::CleanUp()
    Virtual method to eventually reset the specific Event_Phase_Handler.
  */
}

#endif

