#ifndef SHERPA_Single_Events_Colour_Reconnections_H
#define SHERPA_Single_Events_Colour_Reconnections_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/SoftPhysics/Colour_Reconnection_Handler.H"
#include "SHERPA/SoftPhysics/Singlet_Sorter.H"

namespace SHERPA {
  class Colour_Reconnections : public Event_Phase_Handler {
  private :
    bool m_on;
    Colour_Reconnection_Handler * p_reconnectionhandler;
    Singlet_Sorter                m_singlets;
  public :
    Colour_Reconnections(Colour_Reconnection_Handler *);
    ~Colour_Reconnections();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };


  /*!
    \file 
    \brief Contains the class SHERPA::Colour_Reconnections
  */

  /*!
    \class Colour_Reconnections
    \brief The specific Event_Phase_Handler responsible for soft colour reconnections.

    This is the single event phase that preceeds the hadronization stage.
  */
}
#endif

