#ifndef Hard_Processes_H
#define Hard_Processes_H

//#include "XS_Group.H"
#include "XS_Base.H"
#include "Single_Broker.H"
#include "Broker_Group.H"
#include "ISR_Handler.H"
#include "Beam_Handler.H"

#include <string>

namespace APACIC {
  
  class Hard_Processes {
  protected:

    AMEGIC::Process_Base  * two2two;
    ISR::ISR_Handler      * isr;
    BEAM::Beam_Handler    * beam;

    bool ProcessesInit();

  public :
    Hard_Processes(ISR::ISR_Handler *,BEAM::Beam_Handler *,bool &);
    ~Hard_Processes();
    bool OneEvent();
    bool PrepareCalculation();
    bool CalculateCrossSections();

    AMEGIC::Process_Base * proc();
    std::string            Name();

    int                  Nin();
    int                  Nout();
    int               ** Colours();
    AMATOOLS::Vec4D    * Momenta();
    APHYTOOLS::Flavour * Flavs();
  };
  inline AMEGIC::Process_Base * Hard_Processes::proc()    { return two2two; }
  inline std::string            Hard_Processes::Name()    { return std::string("Apacic: Hard_Processes"); }
  inline int                    Hard_Processes::Nin()     { return two2two->Selected()->Nin(); }
  inline int                    Hard_Processes::Nout()    { return two2two->Selected()->Nout(); }
  inline int                 ** Hard_Processes::Colours() { return two2two->Selected()->GetXS()->Colours(); }
  inline AMATOOLS::Vec4D      * Hard_Processes::Momenta() { return two2two->Selected()->Momenta(); }
  inline APHYTOOLS::Flavour   * Hard_Processes::Flavs()   { return two2two->Selected()->Flavs(); }
}
#endif

