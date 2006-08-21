#ifndef Hadron_Decays_H
#define Hadron_Decays_H

#include "Event_Phase_Handler.H"
#include "Hadron_Decay_Handler.H"
#include "Return_Value.H"

namespace ATOOLS {
  class Blob;
}

// #define DEBUG
#ifdef DEBUG
class TH1D;
class TFile;
#endif

namespace SHERPA {
  class Hadron_Decays : public Event_Phase_Handler {
  private :
    HDHandlersMap * p_dechandlers;
#ifdef DEBUG
    TFile* p_file;
    std::map<ATOOLS::kf::code,TH1D> mass_hists;
#endif
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
    \brief Treats blobs containing hadrons to be decayed (which contain status
    needs_hadrondecays)
    */
  /*!
    \fn Hadron_Decays::Treat(ATOOLS::Blob_List *, double &)
    \brief Treats Fragmentation and Hadron Decay Blobs

    This method recurses through the bloblist and when it finds a blob which needs
    hadron decays, it chooses amongst the available decay handlers (e.g. Lund and
    HADRONS, if enabled in Fragmentation.dat). It reiterates until no unstable
    particle is left in a needs_hadrondecays blob.

    While recursing through the unstable particles of one blob, it also puts
    them on a new mass shell. This is done according to the mass distribution
    chosen for this particle in the \c Mass_Handler . If the decay
    handler cannot cope with a chosen mass for the incoming particle, it
    is supposed to deliver a Retry_Method, which will cause this method to dice
    a new mass for this particle and retry it.
    */


}
#endif

