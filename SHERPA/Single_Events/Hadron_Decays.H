#ifndef SHERPA_Single_Events_Hadron_Decays_H
#define SHERPA_Single_Events_Hadron_Decays_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "ATOOLS/Org/Return_Value.H"

namespace ATOOLS {
  class Blob;
  class Particle;
  class Blob_List;
}

namespace SHERPA {
  class Decay_Handler_Base;

  class Hadron_Decays : public Event_Phase_Handler {
  private :
    Decay_Handler_Base * p_dechandler;
    
  public :
    Hadron_Decays(Decay_Handler_Base *);
    ~Hadron_Decays();
    
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
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

    The basic idea is a separation between the creation of the decay blob of
    each daughter particle (CreateDecayBlob), and filling the mother blob with
    kinematics (Treat).
    While the creation step has to make sure
    that everything is known for subsequently dicing the mass of the decaying
    hadron, the filling step does the kinematics in the CMS and then boosts and
    stretches the momenta to the lab CMS and their final masses.
    */
  /*!
    \fn Hadron_Decays::Treat(ATOOLS::Blob_List*)
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

}
#endif

