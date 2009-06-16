#ifndef MODEL_Main_SM_MSSM_EHC_H
#define MODEL_Main_SM_MSSM_EHC_H

#include "MODEL/Main/Model_Base.H"

namespace MODEL {
  class SM_EHC : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
  public :
    SM_EHC(std::string,std::string,bool);
    ~SM_EHC()                   {}
  };

  class MSSM_EHC : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
  public :
    MSSM_EHC(std::string,std::string,bool);
    ~MSSM_EHC()                   {}
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


