#ifndef SHERPA_Single_Events_Beam_Remnants_H
#define SHERPA_Single_Events_Beam_Remnants_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"
#include "SHERPA/SoftPhysics/Fragmentation_Handler.H"

namespace SHERPA {
  class Beam_Remnants : public Event_Phase_Handler {
  private :
    Beam_Remnant_Handler  * p_beamremnanthandler;
    Fragmentation_Handler * p_fragmentationhandler;
  public :
    Beam_Remnants(Beam_Remnant_Handler *);
    ~Beam_Remnants();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };


  /*!
    \file 
    \brief Contains the class SHERPA::Beam_Remnants
  */

  /*!
    \class Beam_Remnants
    \brief The specific Event_Phase_Handler responsible for the translation of partons
           into hadrons.

    This is the single event phase that manages the hadronization stage, i.e. the translation
    of partons into hadrons. For this, it needs some idea on how to treat the remanants of
    eventually coloured incoming beams (like protons) and a fragmentation handler.
  */
  /*!
    \var Fragmentation_Handler * Beam_Remnants::p_fragmentationhandler
    A pointer to the Fragmentation_Handler, filled by the constructor.
  */
  /*!
    \fn Beam_Remnants::Beam_Remnants(Fragmentation_Handler *)
    Copies the pointer to the fragmentation handler and sets the
    name and the type of this event phase. The name is given by 
    "Beam_Remnants : " plus the name of the fragmentation model (e.g. "Lund"), the type
    is "Beam_Remnants".
  */
  /*!
    \fn Beam_Remnants::~Beam_Remnants()
    Nothing to be done in the destructor so far.
  */
  /*!
    \fn bool Beam_Remnants::Treat(ATOOLS::Blob_List *, double &)
    This method basically passes the blob list to the method Hadronize of the fragmentation handler.
  */
  /*!
    \fn void Beam_Remnants::CleanUp()
    Nothing has to be done here.
  */
}
#endif

