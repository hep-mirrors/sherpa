#ifndef Hadronization_h
#define Hadronization_h

#include "Event_Phase_Handler.H"
#include "Beam_Remnant_Handler.H"
#include "Fragmentation_Handler.H"

namespace SHERPA {
  class Hadronization : public Event_Phase_Handler {
  private :
    Beam_Remnant_Handler  * p_beamremnanthandler;
    Fragmentation_Handler * p_fragmentationhandler;
  public :
    Hadronization(Beam_Remnant_Handler *,Fragmentation_Handler *);
    ~Hadronization();
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(std::string) {}
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
    \var Beam_Remnant_Handler * Hadronization::p_beamremnanthandler
    A pointer to the Beam_Remnant_Handler, filled by the constructor.
  */
  /*!
    \var Fragmentation_Handler * Hadronization::p_fragmentationhandler
    A pointer to the Fragmentation_Handler, filled by the constructor.
  */
  /*!
    \fn Hadronization::Hadronization(Beam_Remnant_Handler *,Fragmentation_Handler *)
    Copies the poniters to the beam remnant and the fragmentation handler and sets the
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
    This method basically passes the blob list to the methods FillBeamBlobs and FillBunchBlobs
    of the beam remnant handler and to the method Hadronize of the fragmentation handler.
    These methods fill, if necessary, the blobs for the transition of beams into bunch particles and
    of the bunch particles into partons
  */
  /*!
    \fn void Hadronization::CleanUp()
    Nothing has to be done here.
  */
}
#endif

