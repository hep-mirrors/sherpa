#ifndef SM_MSSM_EHC_h
#define SM_MSSM_EHC_h

#include "Model_Base.H"

namespace MODEL {
  class SM_EHC : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
  public :
    SM_EHC(std::string,std::string,bool);
    ~SM_EHC()                   {}
    bool FillDecay(ATOOLS::Decay_Table *) { return false; }
  };

  class MSSM_EHC : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
  public :
    MSSM_EHC(std::string,std::string,bool);
    ~MSSM_EHC()                   {}
    bool FillDecay(ATOOLS::Decay_Table *) { return false; }
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


