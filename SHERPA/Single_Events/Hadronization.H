#ifndef SHERPA_Single_Events_Hadronization_H
#define SHERPA_Single_Events_Hadronization_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/SoftPhysics/Fragmentation_Handler.H"

namespace SHERPA {
  class Hadronization : public Event_Phase_Handler {
  private :
    Fragmentation_Handler * p_fragmentationhandler;
  public :
    Hadronization(Fragmentation_Handler *);
    ~Hadronization();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };


  /*!
    \file 
    \brief Contains the class SHERPA::Hadronization
  */

  /*!
    \class Hadronization
    \brief The specific Event_Phase_Handler responsible for the translation of partons
           into hadrons.

    This is the single event phase that manages the hadronization stage, i.e. the translation
    of partons into hadrons. For this, it needs some idea on how to treat the remanants of
    eventually coloured incoming beams (like protons) and a fragmentation handler.
  */
  /*!
    \var Fragmentation_Handler * Hadronization::p_fragmentationhandler
    A pointer to the Fragmentation_Handler, filled by the constructor.
  */
  /*!
    \fn Hadronization::Hadronization(Fragmentation_Handler *)
    Copies the pointer to the fragmentation handler and sets the
    name and the type of this event phase. The name is given by 
    "Hadronization : " plus the name of the fragmentation model (e.g. "Lund"), the type
    is "Hadronization".
  */
  /*!
    \fn Hadronization::~Hadronization()
    Nothing to be done in the destructor so far.
  */
  /*!
    \fn bool Hadronization::Treat(ATOOLS::Blob_List *, double &)
    This method basically passes the blob list to the method Hadronize of the fragmentation handler.
  */
  /*!
    \fn void Hadronization::CleanUp()
    Nothing has to be done here.
  */
}
#endif

