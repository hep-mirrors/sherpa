#ifndef MODEL_Main_SM_U1_B_H
#define MODEL_Main_SM_U1_B_H

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"

namespace MODEL {
  class Standard_Model;
  class SM_U1_B : public Model_Base {
  private :
    Standard_Model *p_sm;
    void ParticleInit();
    void FixMix();
    void FillSpectrum(PDF::ISR_Handler *const isr);
    bool SanityChecks();
  public :
    SM_U1_B(std::string,std::string,bool);
    ~SM_U1_B();
    bool ModelInit(PDF::ISR_Handler *const isr);
  };
  /*!
    \class SM_U1_B
    \brief contains the SM_U1_B (no surprise)
  */
  /*!
    \var std::string SM_U1_B::m_generator
    The spectrum generator to be used: None.
  */
  /*!
    \fn SM_U1_B::SM_U1_B(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, 
    MODEL::SM.  Its parameters (i.e. constants, functions, matrices) 
    are then copied into the newly initialized pointers to the 
    searchable maps. 
  */
  /*!
    \fn void SM_U1_B::FillSpectrum()
    This method is called by the constructor and selects the model framework. 
    There are two parameters (in addition to the SM ones):
    - \f$M_{Z'}f$, the mass of the Z'
    - \f$g'_1\f$, the gauge coupling of the additional \f$U(1)\f$, which, 
    however only couples to baryon number (i.e. democratically to quarks).
  */
}

#endif


