#ifndef Shower_Initialization_h
#define Shower_Initialization_h

#include "Model_Base.H"
#include "ISR_Handler.H"
#include "Hard_Interface.H"

#include <string>
#include <iostream> 

namespace SHERPA {
  class Shower_Initialization {
  private :
    std::string                  m_dir,m_file;
    int                          m_maxjetnumber;

    AORGTOOLS::Data_Read       * p_dataread;

    MODEL::Model_Base          * p_model;
    ISR::ISR_Handler           * p_isr;
    APACIC::Hard_Interface     * p_apacic;
  private :
    void ReadInFile();
    bool InitializeApacic(bool,bool);
  public :
    Shower_Initialization(std::string,std::string,MODEL::Model_Base *,
			  ISR::ISR_Handler *,int);
    ~Shower_Initialization();

    void Output();
  };
}


#endif
