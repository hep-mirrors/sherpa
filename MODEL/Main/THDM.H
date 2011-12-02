#ifndef MODEL_Main_THDM_H
#define MODEL_Main_THDM_H

#include "MODEL/Main/Model_Base.H"

namespace MODEL {
  class Standard_Model;
  class THDM : public Model_Base {
  private :
    Standard_Model *p_sm;
    void ParticleInit();
    void FillSpectrum(const PDF::ISR_Handler_Map& isr);
  public :
    THDM(std::string,std::string,bool);
    ~THDM();
    bool ModelInit(const PDF::ISR_Handler_Map& isr);
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


