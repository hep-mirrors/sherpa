#ifndef MODEL_Main_SM_AGC_H
#define MODEL_Main_SM_AGC_H

#include "MODEL/Main/Model_Base.H"

namespace MODEL {
  class SM_AGC : public Model_Base {
  private :
    void ParticleInit() {}
    void FillSpectrum();
  public :
    SM_AGC(std::string,std::string,bool);
    ~SM_AGC()                   {}
  };
  /*!
    \class SM_AGC
    \brief contains the Standard_Model extended by anomalous triple & quartic
    gauge couplings (see hep-ph/0001065 and Nucl. Phys. B282 (1987) 253-307)
    
  */
  /*!
    \fn void SM_AGC::FillSpectrum()
    This method is responsible for setting all the parameters of the model that
    have to be specified in Model.dat
  */
}

#endif


