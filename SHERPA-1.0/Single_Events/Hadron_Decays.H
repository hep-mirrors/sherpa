#ifndef Hadron_Decays_H
#define Hadron_Decays_H

#include "Event_Phase_Handler.H"
#include "Hadron_Decay_Handler.H"

namespace SHERPA {
  class Hadron_Decays : public Event_Phase_Handler {
  private :
    HDHandlersMap * p_dechandlers;
  public :
    Hadron_Decays(HDHandlersMap *);
    ~Hadron_Decays();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp();
    void                       Finish(const std::string &);
  };

  /*!
    \file Hadron_Decays.H
    \brief Declares the class SHERPA::Hadron_Decays
    */
  /*!
    \class Hadron_Decays
    \brief Treats blobs containing hadrons to be decayed
    */
  /*!
    \fn Hadron_Decays::Treat(ATOOLS::Blob_List *, double &)
    \brief Treats Fragmentation and Hadron Decay Blobs

    This method calls the hadron decay handler to treat a blob that contains
    unstable hadron: Hadron_Decay_Handler::FillHadronDecayBlobs.
    It also decides if it has to call the Pythia or the HADRONS hadron decay handler
    
    Basically, it goes through the blob list, picks unstable hadrons and makes them decay
    until all unstable hadrons are developed.
    */


}
#endif

