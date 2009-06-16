#ifndef CSSHOWER_Sudakov_CS_IntegratedSudakov_H
#define CSSHOWER_Sudakov_CS_IntegratedSudakov_H

#include "CSSHOWER++/Sudakov/CS_SingleSudakov_Base.H"
#include <list>

namespace MODEL { class Running_AlphaS; }
namespace PDF   { 
  class ISR_Handler; 
  class PDF_Base; 
}

namespace CSSHOWER {
  
  class CS_IntegratedSudakov {
  
  private:

    std::list<CS_SingleSudakov_Base *> m_suds;
    std::list<CS_SingleSudakov_Base *>::iterator m_ssud;
    
    int                     m_mode;
    MODEL::Running_AlphaS * p_runas;
    PDF::ISR_Handler      * p_isr;
    PDF::PDF_Base         * p_pdf[2];
    
  public:

    // constructor 
    CS_IntegratedSudakov(int mode,MODEL::Running_AlphaS *as,PDF::ISR_Handler * isr); 
    // destructor
    ~CS_IntegratedSudakov();
    
    //methods
    void   Add(CS_SingleSudakov_Base *);
    double Delta(const cstp::code dtype,const ATOOLS::Flavour & fl,
		 const double Q2, const double mk2,const double kt2prod,
		 const double kt2dec);
  };
}

#endif
