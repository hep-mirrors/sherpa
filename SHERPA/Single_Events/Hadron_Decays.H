#ifndef SHERPA_Single_Events_Hadron_Decays_H
#define SHERPA_Single_Events_Hadron_Decays_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

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
  class Soft_Photon_Handler;

  class Hadron_Decays : public Event_Phase_Handler {
  private :
#ifdef DEBUG__Hadrons
#ifdef USING__ROOT
    TFile* p_file;
    std::map<kf_code,TH1D*> mass_hists;
#endif
#endif
    HDHandlersMap     * p_dechandlers;
    Soft_Photon_Handler * p_softphotons;
    ATOOLS::Blob_List * p_bloblist;
    
    static int          m_mass_smearing;
    
    bool Treat(ATOOLS::Blob * blob);
    bool CreateDecayBlob(ATOOLS::Particle* inpart);
    bool SetMasses(ATOOLS::Blob*);
    bool BoostAndStretch(ATOOLS::Blob* blob, const ATOOLS::Vec4D& labmom);
    Hadron_Decay_Handler* ChooseDecayHandler(ATOOLS::Particle* part);
    bool RejectExclusiveChannelsFromFragmentation(ATOOLS::Blob*);
    bool PerformMixing(ATOOLS::Particle* inpart, ATOOLS::Blob_List* bloblist);
    bool AttachExtraQED(ATOOLS::Blob* blob);
    void SetPosition(ATOOLS::Blob* blob);
    
  public :
    Hadron_Decays(HDHandlersMap *, Soft_Photon_Handler *);
    ~Hadron_Decays();
    
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp();
    void                       Finish(const std::string &);

    static void SetMassSmearing(int mode) { m_mass_smearing=mode; }
    static int MassSmearing()             { return m_mass_smearing; }
  };

  /*!
    \file Hadron_Decays.H
    \brief Declares the class SHERPA::Hadron_Decays
    */
  /*!
    \class Hadron_Decays
    \brief Treats blobs containing hadrons to be decayed (which contain status
    needs_hadrondecays)

    The basic idea is a separation between the creation of the decay blob of
    each daughter particle (CreateDecayBlob), and filling the mother blob with
    kinematics (Treat).
    While the creation step has to make sure
    that everything is known for subsequently dicing the mass of the decaying
    hadron, the filling step does the kinematics in the CMS and then boosts and
    stretches the momenta to the lab CMS and their final masses.
    */
  /*!
    \var int Hadron_Decays::m_mass_smearing
    \brief Switch to turn mass smearing on and off.

    Can be set by the "MASS_SMEARING = [0,1]" variable in Fragmentation.dat
  */
  /*!
    \fn Hadron_Decays::Treat(ATOOLS::Blob_List *, double &)
    \brief Treats all blobs in a blob list which need hadron decays.

    This method recurses through the bloblist and when it finds a blob which needs
    hadron decays, it decays all its daughters. Care is taken for the following
    special features:
    - The daughters are treated in random order to avoid biases in spin correlations
    - If a fragmentation blob is encountered, it is checked whether it stems from a
      partonic hadron decay. If it does, and it produced an exclusive decay channel
      of that decaying hadron, it is rejected and redone, as to not spoil that
      hadrons branching ratio.
    - all outparticles of the incoming blob are set off-shell according to their
      decay channels.
    */
  /*!
    \fn Hadron_Decays::Treat(ATOOLS::Blob*)
    \brief Treats the given blob by adding decay products kinematics, and
    boosting and stretching the blob such that it finally is back in the lab
    frame and the daughter particles have off-shell masses.
    For that it first has to do the on-shell kinematics, then create the
    daughter decay blobs to know how to smear the daughter masses, and finally
    boost and stretch everything back.
    */
  /*!
    \fn Hadron_Decays::CreateDecayBlob(ATOOLS::Particle* inpart)
    \brief Chooses the appropriate decay handler for <var>inpart</var>
    to create the decay blob for inpart. The decay blob is only a stub
    so far without kinematics and without outparticles.
    */
  /*!
    \fn Hadron_Decays::SetMasses(ATOOLS::Blob*)
    \brief Uses the appropriate decay handler for each outparticle to
    generate its destination mass, potentially using the information attached
    in the CreateDecayBlob method. The generated masses are stored in the
    FinalMass() property of each particle.
    */
  /*!
    \fn Hadron_Decays::BoostAndStretch(ATOOLS::Blob*)
    \brief  This method does multiple things with a CMS & on-shell blob:
    1. set inparticle momentum off-shell to the labmom mass (== FinalMass())
    2. set the outparticle momenta off-shell to their FinalMass()
    3. boost the blob into the lab as defined by labmom
    */
  /*!
    \fn Hadron_Decay_Handler* ChooseDecayHandler(ATOOLS::Particle* part)
    \brief Chooses a hadron decay handler which can deal with the
    particle according to their priorities.
    */
  /*!
    \fn bool RejectExclusiveChannelsFromFragmentation(ATOOLS::Blob*)
    \brief If a fragmentation blob is encountered, it is checked whether it stems from a
    partonic hadron decay. If it does, and it produced an exclusive decay channel
    of that decaying hadron, it is rejected and redone, as to not spoil that
    hadron's branching ratio.
    */
  /*!
    \fn PerformMixing(ATOOLS::Particle* inpart, ATOOLS::Blob_List* bloblist)
    \brief Uses the appropriate hadron decay handler for inpart to add a
    mixing blob for it in case that is activated and generated.
    */
  /*!
    \fn void SetPosition(ATOOLS::Blob* blob)
    \brief Calculates the position of the blob in the lab frame according
    to its inparticle's decay time and momentum.
    */

}
#endif

