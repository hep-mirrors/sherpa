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
  };
  /*!
    \class Hadronization
    \brief
  */
  /*!
    \var Beam_Remnant_Handler * Hadronization::p_beamremnanthandler
  */
  /*!
    \var Fragmentation_Handler * Hadronization::p_fragmentationhandler
  */
  /*!
    \fn Hadronization::Hadronization(Beam_Remnant_Handler *,Fragmentation_Handler *)
  */
  /*!
    \fn Hadronization::~Hadronization()
  */
  /*!
    \fn bool Hadronization::Treat(ATOOLS::Blob_List *, double &)
  */
  /*!
    \fn void Hadronization::CleanUp()
    Nothing has to be done here.
  */
}
#endif

