#ifndef Hadron_Decay_Handler_h
#define Hadron_Decay_Handler_h


#include <set> 
#include <iterator> 
#include <iostream> 

#include "Lund_Interface.H"
#include "Hadron_Decays.H"
#include "Particle_List.H"
#include "Particle.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Hadron_Decay_Handler {
  protected :
    std::string              m_dir,m_file,m_decmodel;
    int                      m_mode;
    Lund_Interface         * p_lund; 
    HADRONS::Hadron_Decays * p_hadrons;

    bool ReconstructLundHadronDecays();

  public :

    Hadron_Decay_Handler(std::string,std::string,HADRONS::Hadron_Decays *);
    Hadron_Decay_Handler(std::string,std::string,Lund_Interface *);

    ~Hadron_Decay_Handler();

    bool FillHadronDecayBlobs(ATOOLS::Blob_List *,ATOOLS::Particle_List *);

  };

  /*!
    \class Hadron_Decay_Handler
    \brief So far not being used.
    
    This class so far is not being used. Work is ongoing.
  */
}
#endif

