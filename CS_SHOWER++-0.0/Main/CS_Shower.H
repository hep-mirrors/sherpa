#ifndef CS_Shower_h
#define CS_Shower_h

#include "Shower.H"
#include "Singlet.H"
#include "Blob_List.H"
#include "Blob.H"

namespace PDF    { class ISR_Handler;}

namespace CS_SHOWER {
  class CS_Shower {
  private:
    bool m_isron,m_fsron;
    All_Singlets    * p_allsinglets;
    Shower            m_shower;
    ATOOLS::Vec4D   * p_momenta;
    ATOOLS::Flavour * p_flavours;
    ATOOLS::Poincare  m_boost, m_zrot;
    bool InitBreitFrame(ATOOLS::Blob *);
    void BoostBack(ATOOLS::Vec4D&);
    void ExtractPartons(ATOOLS::Blob *,ATOOLS::Blob *);
  public:
    CS_Shower(PDF::ISR_Handler * isr, bool isron,bool fsron);
    ~CS_Shower();
    bool PerformShowers(All_Singlets * =NULL);
    bool ExtractPartons(ATOOLS::Blob_List *);
    void PrepareAllSinglets();
    bool ISROn() {return m_isron;}            
    bool FSROn() {return m_fsron;}            
    
    inline All_Singlets * GetAllSinglets() { return p_allsinglets; }
  };
}

#endif
