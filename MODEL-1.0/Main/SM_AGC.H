#ifndef SM_AGC_h
#define SM_AGC_h

#include "Model_Base.H"

namespace MODEL {
  class SM_AGC : public Model_Base {
  private :
    void FillSpectrum();
  public :
    SM_AGC(std::string,std::string);
    ~SM_AGC()                   {}
    bool FillDecay(ATOOLS::Decay_Table *) { return false; }
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


