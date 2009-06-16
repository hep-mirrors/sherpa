#ifndef SHERPA_SoftPhysics_Hadron_Decay_Handler_H
#define SHERPA_SoftPhysics_Hadron_Decay_Handler_H

#include <set>
#include <iterator>

#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/CXXFLAGS.H"


namespace HADRONS { class Hadrons; }
namespace ATOOLS { class Blob; class Particle; class Blob_List; class Cluster_Amplitude; }

namespace SHERPA {
#ifdef USING__PYTHIA
  class Lund_Interface;
#endif

  class Hadron_Decay_Handler: public ATOOLS::Mass_Selector {
  protected :
    std::string                   m_decmodel;
    int                           m_mode;
    HADRONS::Hadrons            * p_hadrons;
#ifdef USING__PYTHIA
    Lund_Interface              * p_lund; 
#endif
    std::set<kf_code>  * p_cans;
    ATOOLS::Cluster_Amplitude   * p_ampl;

  public :
    Hadron_Decay_Handler(HADRONS::Hadrons *);
#ifdef USING__PYTHIA
    Hadron_Decay_Handler(Lund_Interface *);
#endif
    ~Hadron_Decay_Handler();
    
    bool CreateDecayBlob(ATOOLS::Blob*);
    bool FillDecayBlob(ATOOLS::Blob *, const ATOOLS::Vec4D& labmom);
    
    bool CanDealWith(kf_code);

    ATOOLS::Cluster_Amplitude *ClusterConfiguration(ATOOLS::Blob *const bl);

    HADRONS::Hadrons * GetHadrons() { return p_hadrons; }
#ifdef USING__PYTHIA
    Lund_Interface   * GetLund()    { return p_lund; }
#endif
    std::string        Name()       { return m_decmodel; }
    std::set<kf_code>  * GetCans() { return p_cans; }
    bool               DiceMass(ATOOLS::Particle*, double min, double max);
    void               SetSignalProcessBlob(ATOOLS::Blob* spblob);
    bool               PerformMixing(ATOOLS::Particle* inpart, ATOOLS::Blob_List* bloblist);
    void               CleanUp();
    bool               IsExclusiveDecaychannel(ATOOLS::Blob* blob, ATOOLS::FlavourSet results);

    double Mass(const ATOOLS::Flavour &fl) const;

  };

  typedef std::map<std::string,Hadron_Decay_Handler *> HDHandlersMap;
  typedef HDHandlersMap::iterator                      HDHandlersIter;
}

  /*!
    \file Hadron_Decay_Handler.H
    \brief Declares the class SHERPA::Hadron_Decay_Handler
    */
  /*!
    \class Hadron_Decay_Handler
    \brief A wrapper class around all implemented hadron decay modules.

    Knows all different hadron decay handlers and can tell whether a specific
    one can deal with an unstable particle (CanDealWith()). It can then pass
    the to be decayed blob to it (CreateDecayBlob() and FillDecayBlob()).
    */
  /*!
    \fn CreateDecayBlob(ATOOLS::Blob *)
    \brief Passes the to be created blob to the hadron decay handler
    incorporated.
    The handler then has the chance to prepare the blob and attach any
    information necessary for its inparticle's mass dicing to it.
    */
  /*!
    \fn FillDecayBlob(ATOOLS::Blob *)
    \brief Passes the to be filled blob to the hadron decay handler
    incorporated. It calls the decay method of its hadron decay handler
    for the given blob. The blob is in the incoming particle's rest frame
    all the time.
    */

#endif
