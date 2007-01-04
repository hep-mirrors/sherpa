#ifndef Hadron_Decays_H
#define Hadron_Decays_H

#include "Event_Phase_Handler.H"
#include "Hadron_Decay_Handler.H"
#include "Return_Value.H"

namespace ATOOLS {
  class Blob;
  class Particle;
  class Blob_List;
  class Amplitude_Tensor;
}

#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
#include "TH1D.h"
#include "TFile.h"
#endif
#endif

namespace SHERPA {
  class Hadron_Decays : public Event_Phase_Handler {
  private :
    HDHandlersMap            * p_dechandlers;
    ATOOLS::Blob_List        * p_bloblist;
    ATOOLS::Particle_Vector    m_daughters;
    std::vector<ATOOLS::Vec4D> m_saved_momenta;
    ATOOLS::Amplitude_Tensor*  p_saved_amplitudes;
    double                     m_max_mass;
    bool                       m_mass_smearing;
#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
    TFile* p_file;
    std::map<ATOOLS::kf::code,TH1D*> mass_hists;
#endif
#endif

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob * blob);
    ATOOLS::Return_Value::code Treat(int i);
    bool PrepareMassSmearing(ATOOLS::Blob* blob);
    void ResetMotherAmplitudes(ATOOLS::Blob* blob);
    Hadron_Decay_Handler* ChooseDecayHandler(ATOOLS::Particle* part);
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
    \var double Hadron_Decays::m_mass_smearing
    \brief Switch to turn mass smearing on and off.

    Can be set by the "MASS_SMEARING = [0,1]" variable in Fragmentation.dat
  */
  /*!
    \var double Hadron_Decays::m_max_mass
    \brief Temporary storage for total incoming mass in blob.
  */
  /*!
    \fn Hadron_Decays::Treat(ATOOLS::Blob_List *, double &)
    \brief Treats all blobs in a blob list which need hadron decays.

    This method recurses through the bloblist and when it finds a blob which needs
    hadron decays, it treats it with the Treat(Blob*) method.
    */
  /*!
    \fn Hadron_Decays::Treat(ATOOLS::Blob*)
    \brief Treats the given blob by decaying its outgoing hadrons.

    It uses the Treat(int i) method to decay all outgoing particles. If this
    method returns Failure for one particle, it starts over from the beginning.
    If 10 retries didn't help, it returns Retry_Event, otherwise Success.
    */
  /*!
    \fn Hadron_Decays::Treat(int i)
    \brief Decays the daughter particle #i.

    It chooses amongst the available decay handlers (e.g. Lund and
    HADRONS, if enabled in Fragmentation.dat).

    It also puts the particle on a new mass shell. This is done according to
    the mass distribution chosen for this particle in the \c Mass_Handler .
    If the decay handler cannot cope with a chosen mass for the incoming
    particle, it is supposed to deliver a Retry_Method, which will cause this
    method to dice a new mass for this particle and retry it. If 100 retries
    didn't help, it returns Failure.
    */


}
#endif

