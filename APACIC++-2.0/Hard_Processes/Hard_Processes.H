#ifndef Hard_Processes_H
#define Hard_Processes_H

//#include "XS_Group.H"
#include "XS_Base.H"
#include "All_Processes.H"
#include "ISR_Handler.H"
#include "Beam_Spectra_Handler.H"

#include <string>

namespace APACIC {
  
  class Hard_Processes {
  protected:

    EXTRAXS::XS_Base         * two2two;

    AMEGIC::All_Processes    * broker;
    APHYTOOLS::Selector_Data * seldata;
    ISR::ISR_Handler         * isr;
    BEAM::Beam_Spectra_Handler       * beam;

    bool ProcessesInit();

  public :
    Hard_Processes(APHYTOOLS::Selector_Data *,ISR::ISR_Handler *,BEAM::Beam_Spectra_Handler *,bool &);
    ~Hard_Processes();
    bool OneEvent();
    bool PrepareCalculation();
    bool CalculateCrossSections();

    //    AMEGIC::Process_Base * GetBroker();
    std::string            Name();

    int                  Nin();
    int                  Nout();
    int               ** Colours();
    AMATOOLS::Vec4D    * Momenta();
    APHYTOOLS::Flavour * Flavs();
  };
  //  inline AMEGIC::Process_Base * Hard_Processes::GetBroker() { return broker; }
  inline std::string            Hard_Processes::Name()      { return std::string("Apacic: Hard_Processes"); }
  inline int                    Hard_Processes::Nin()       { return two2two->Nin();      }
  inline int                    Hard_Processes::Nout()      { return two2two->Nout();     }
  inline int                 ** Hard_Processes::Colours()   { return two2two->Colours();  }
  inline AMATOOLS::Vec4D      * Hard_Processes::Momenta()   { return broker->Selected()->Momenta(); }
  inline APHYTOOLS::Flavour   * Hard_Processes::Flavs()     { return two2two->Flavs();    }
}
#endif

