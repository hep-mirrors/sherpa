#ifndef SHERPA_PerturbativePhysics_Shower_Handler_H
#define SHERPA_PerturbativePhysics_Shower_Handler_H

#include "ATOOLS/Phys/Blob_List.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/Shower_Base.H"

namespace MODEL    { class Model_Base; }
namespace REMNANTS { class Remnants_Handler; }

namespace SHERPA {

  class Shower_Handler {
  private:
    PDF::Shower_Base *p_shower;
    PDF::ISR_Handler *p_isr;       
    std::string m_name;
  public:
    //constructor
    Shower_Handler(MODEL::Model_Base *const,
                   PDF::ISR_Handler *const,
                   const int isrtype);

    // destructor
    ~Shower_Handler();

    // member functions
    void  SetRemnants(REMNANTS::Remnant_Handler*);
    void  FillBlobs(ATOOLS::Blob_List *const bl);
    void  FillDecayBlobs(ATOOLS::Blob_List *const bl);
    void  CleanUp();

    // inline functions
    inline PDF::Shower_Base  * GetShower()     { return p_shower; }
    inline PDF::ISR_Handler  * GetISRHandler() { return p_isr;    }
    inline std::string ShowerGenerator()       { return m_name; }
    inline bool On()                           { return p_shower!=NULL; }

  };
  typedef std::map<PDF::isr::id, Shower_Handler*> Shower_Handler_Map;

}

#endif
