#ifndef Output_Handler_h
#define Output_Handler_h

#include "HepMC_Interface.H"
#include "Blob_List.H"

namespace SHERPA {
  class Output_Handler {
  private :
    HepMC_Interface * p_hepmc;
    HepMC::GenEvent * p_event;
    bool              m_active;
    int               m_type;
  public :
    Output_Handler(int type=1);
    ~Output_Handler();
    void OutputToFormat(APHYTOOLS::Blob_List *);
    bool Active() { return m_active; }
  };
}
#endif
