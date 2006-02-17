#ifndef Hadron_Decay_Handler_h
#define Hadron_Decay_Handler_h


#include <set> 
#include <iterator> 
#include <iostream> 

#include "Lund_Interface.H"
#include "Particle_List.H"
#include "Particle.H"
#include "Blob_List.H"
#include "Blob.H"
#include "CXXFLAGS.H"

#ifdef USING__Hadrons
#include "Hadrons.H"
#endif
#include "Spin_Density_Matrix.H"

namespace SHERPA {
  class Hadron_Decay_Handler {
  protected :
    std::string              m_decmodel;
    int                      m_mode;
#ifdef USING__Hadrons
    HADRONS::Hadrons       * p_hadrons;
#endif
    Lund_Interface         * p_lund; 

    void SwitchOfLundDecays();

  public :
#ifdef USING__Hadrons
    Hadron_Decay_Handler(HADRONS::Hadrons *);
    HADRONS::Hadrons * GetHadrons() { return p_hadrons; }
#endif
    Hadron_Decay_Handler(Lund_Interface *);

    ~Hadron_Decay_Handler();
	void DeletePointers(); 

    void PrepareDecays(ATOOLS::Blob *);
    bool FillHadronDecayBlobs(ATOOLS::Particle *,ATOOLS::Blob_List *,
                  ATOOLS::Spin_Density_Matrix * sigma = NULL, 
                  ATOOLS::Spin_Density_Matrix * decmatr = NULL,
			      ATOOLS::Particle_List *part_list=NULL );
    Lund_Interface * GetLund() { return p_lund; }
    void EraseTreated(std::set<int> *);
  };

  typedef std::map<std::string,Hadron_Decay_Handler *> HDHandlersMap;
  typedef HDHandlersMap::iterator                      HDHandlersIter;
 
  /*!
    \file Hadron_Decay_Handler.H
    \brief Declares the class SHERPA::Hadron_Decay_Handler
    */
  /*!
    \class Hadron_Decay_Handler
    \brief Redirects a decay to the corresponding module (Pythia, Hadrons, ...)
    
    This class knows the implemented hadron decay modules and tells them to make
    an unstable hadron decay if the Hadron_Decays::Treat method tells them to do so.
  */
  /*!
    \fn Hadron_Decay_Handler::FillHadronDecayBlobs(ATOOLS::Particle *,ATOOLS::Blob_List *, ATOOLS::Spin_Density_Matrix * sigma = NULL, ATOOLS::Spin_Density_Matrix * decmatr = NULL, ATOOLS::Particle_List *part_list=NULL )
    \brief Calls the method that actually performs the decay (chain)

    This method calls either the PerformDecay method of Pythia or the Hadrons module.
    */
  /*!
    \fn Hadron_Decay_Handler::EraseTreated(std::set<int> *)
    This method goes through all known hadron decay handlers (HADRONS) and check with which
    unstable hadrons they can cope with. These are erased from the list, so that at the end
    the Initialization_Handler::InitializeTheHadronDecays knows which particles Sherpa can make decay
    and which it cannot. Furthermore, it tells the Spin Correlation Tensor which particles needs to
    be considered (only the Sherpa hadron decay handler, HADRONS, has an implemented spin correlation algorithm).
    */
  /*!
    \fn Hadron_Decay_Handler::SwitchOfLundDecays()
    This methods tells Pythia to set those hadrons stable which HADRONS can cope with.
    */

}
#endif

