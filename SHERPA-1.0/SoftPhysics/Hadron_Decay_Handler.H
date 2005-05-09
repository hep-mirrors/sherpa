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

namespace HADRONS { class Hadrons; }

namespace SHERPA {
  class Hadron_Decay_Handler {
  protected :
    std::string              m_decmodel;
    int                      m_mode;
    Lund_Interface         * p_lund; 
    HADRONS::Hadrons       * p_hadrons;

	void SwitchOfLundDecays();

  public :

    Hadron_Decay_Handler(HADRONS::Hadrons *);
    Hadron_Decay_Handler(Lund_Interface *);

    ~Hadron_Decay_Handler();
	void DeletePointers(); 

    bool FillHadronDecayBlobs(ATOOLS::Particle *,ATOOLS::Blob_List *,
							  ATOOLS::Particle_List *part_list=NULL);
	HADRONS::Hadrons * GetHadrons() { return p_hadrons; }
	Lund_Interface * GetLund() { return p_lund; }
	void EraseTreated(std::set<int> *);
  };

  typedef std::map<std::string,Hadron_Decay_Handler *> HDHandlersMap;
  typedef HDHandlersMap::iterator                      HDHandlersIter;
 
  /*!
    \class Hadron_Decay_Handler
    \brief So far not being used.
    
    This class so far is not being used. Work is ongoing.
  */
}
#endif

