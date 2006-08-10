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
    std::string                   m_decmodel;
    int                           m_mode;
#ifdef USING__Hadrons
    HADRONS::Hadrons            * p_hadrons;
#endif
    Lund_Interface              * p_lund; 
    std::set<ATOOLS::kf::code>  * p_cans;

    void SwitchOfLundDecays();
  public :
#ifdef USING__Hadrons
    Hadron_Decay_Handler(HADRONS::Hadrons *);
#endif
    Hadron_Decay_Handler(Lund_Interface *);

    ~Hadron_Decay_Handler();

    bool FillHadronDecayBlob(ATOOLS::Particle *,ATOOLS::Blob_List *);
    void EraseTreated(std::set<int> *);
    bool CanDealWith(ATOOLS::kf::code);

#ifdef USING__Hadrons
    HADRONS::Hadrons * GetHadrons() { return p_hadrons; }
#endif
    Lund_Interface   * GetLund()    { return p_lund; }
    std::string        Name()       { return m_decmodel; }
  };

  typedef std::map<std::string,Hadron_Decay_Handler *> HDHandlersMap;
  typedef HDHandlersMap::iterator                      HDHandlersIter;
}
#endif

