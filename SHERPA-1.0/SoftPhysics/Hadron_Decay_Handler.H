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

namespace SHERPA {
  class Hadron_Decay_Handler {
  protected :

    std::string    m_dir,m_file;
    Lund_Interface * p_lund; 

    bool ReconstructLundHadronDecays();

  public :

    Hadron_Decay_Handler(std::string,std::string);
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

