#ifndef Shower_Handler_h
#define Shower_Handler_h 

#include "Model_Base.H"
#include "ISR_Handler.H"
#include "Hard_Interface.H"
#include "Data_Read.H"

#include <string>
#include <iostream> 

namespace APACIC { class Tree; }

namespace SHERPA {
  class Shower_Handler {
  private :
    std::string                  m_dir,m_file;
    std::string                  m_showergenerator;
    int                          m_maxjetnumber;
    bool                         m_isrshowerswitch,m_fsrshowerswitch;

    AORGTOOLS::Data_Read       * p_dataread;
    MODEL::Model_Base          * p_model;
    ISR::ISR_Handler           * p_isr;
    APACIC::Hard_Interface     * p_apacic;
  public :
    Shower_Handler(std::string,std::string,MODEL::Model_Base *,
		   ISR::ISR_Handler *,int);
    ~Shower_Handler();
    int             PerformShowers(bool);
    void            FillBlobs(APHYTOOLS::Blob_List *);
    void            CleanUp();
    std::string     ShowerGenerator() { return m_showergenerator; }
    bool            ISROn()           { return m_isrshowerswitch; }
    bool            FSROn()           { return m_fsrshowerswitch; }
    int             MaxJetNumber()    { return m_maxjetnumber;    }
    APACIC::Tree  * GetFinTree(); 
    APACIC::Tree ** GetIniTrees();
  };
}


#endif
