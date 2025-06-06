#ifndef SHERPA_Single_Events_Nuclear_Formation_H
#define SHERPA_Single_Events_Nuclear_Formation_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/SoftPhysics/Colour_Reconnection_Handler.H"
#include "ATOOLS/Phys/Fragmentation_Base.H"
#include "SHERPA/SoftPhysics/Singlet_Sorter.H"

namespace SHERPA {
  class Nuclear_Formation : public Event_Phase_Handler {
  private :
    bool m_on;
    Colour_Reconnection_Handler * p_reconnectionhandler;
    ATOOLS::Fragmentation_Base  * p_fragmentationhandler;
    Singlet_Sorter                m_singlets;
  public :
    Nuclear_Formation(Colour_Reconnection_Handler *,
		  ATOOLS::Fragmentation_Base *);
    ~Nuclear_Formation();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };


  /*!
    \file 
    \brief Contains the class SHERPA::Nuclear_Formation
  */

  /*!
    \class Nuclear_Formation
    \brief The specific Event_Phase_Handler responsible for nucleus formation

    This is the single event phase that manages the hadronization stage, i.e. the translation
    of nucleons into nucleus. For this, it needs some idea on how to treat the remanants of
    eventually coloured incoming beams (like protons) and a fragmentation handler.
  */
  /*!
    \var Fragmentation_Handler * Hadronization::p_fragmentationhandler
    A pointer to the Fragmentation_Handler, filled by the constructor.
  */
  /*!
    \fn Nuclear_Formation::Nuclear_Formation(Fragmentation_Handler *)
    Copies the pointer to the fragmentation handler and sets the
    name and the type of this event phase. The name is given by 
    "Nuclear_Formation : " plus the name of the fragmentation model, the type is "Nuclear_Formation".
  */
  /*!
    \fn Nuclear_Formation::~Nuclear_Formation()
    Nothing to be done in the destructor so far.
  */
  /*!
    \fn bool Nuclear_Formation::Treat(ATOOLS::Blob_List*)
    This method basically passes the blob list to the method Hadronize of the fragmentation handler.
  */
  /*!
    \fn void Nuclear_Formation::CleanUp()
    Nothing has to be done here.
  */
}
#endif

