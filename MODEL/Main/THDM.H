#ifndef THDM_h
#define THDM_h

#include "Model_Base.H"

namespace MODEL {
  class THDM : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
  public :
    THDM(std::string,std::string,bool);
    ~THDM();
    bool FillDecay(ATOOLS::Decay_Table *) { return 0; }
  };
  /*!
    \class THDM
    \brief contains the THDM (no surprise)
  */
  /*!
    \fn THDM::THDM(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, MODEL::SM. 
    Its parameters (i.e. constants, functions, matrices) are then copied into the newly 
    initialized pointers to the searchable maps. 
  */
  /*!
    \fn void THDM::FillSpectrum()
    This method is called by the constructor and selects the model framework.
  */
}

#endif


