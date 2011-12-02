#ifndef MODEL_Main_SM_MSSM_EHC_H
#define MODEL_Main_SM_MSSM_EHC_H

#include "MODEL/Main/Model_Base.H"

namespace MODEL {
  class Standard_Model;
  class SM_EHC : public Model_Base {
  private :
    Standard_Model *p_sm;
    void ParticleInit();
    void FillSpectrum(const PDF::ISR_Handler_Map& isr);
  public :
    SM_EHC(std::string,std::string,bool);
    ~SM_EHC()                   {}
    bool ModelInit(const PDF::ISR_Handler_Map& isr);
  };

  class MSSM;
  class MSSM_EHC : public Model_Base {
  private :
    MSSM *p_mssm;
    void ParticleInit();
    void FillSpectrum(const PDF::ISR_Handler_Map& isr);
  public :
    MSSM_EHC(std::string,std::string,bool);
    ~MSSM_EHC()                   {}
    bool ModelInit(const PDF::ISR_Handler_Map& isr);
  };

  /*!
    \class SM_EHC
    \brief contains the Standard_Model extended by a loop induced Gluon-Gluon-Higgs coupling
    
  */
  /*!
    \fn void SM_EHC::FillSpectrum()
    This method is responsible for setting all the parameters of the model.
  */

  /*!
    \class MSSM_EHC
    \brief contains the MSSM extended by a loop induced Gluon-Gluon-Higgs couplings
    
  */
  /*!
    \fn void MSSM::FillSpectrum()
    This method is responsible for setting all the parameters of the model.
  */
}

#endif


