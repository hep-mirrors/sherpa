#ifndef MODEL_Main_SM_Phantom_U1_H
#define MODEL_Main_SM_Phantom_U1_H

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"

namespace MODEL {
  class Standard_Model;
  class SM_Phantom_U1 : public Model_Base {
  private :
    Standard_Model *p_sm;
    void ParticleInit();
    void FillSpectrum();
    bool SanityChecks();
  public :
    SM_Phantom_U1(std::string,std::string,bool);
    ~SM_Phantom_U1();
    bool ModelInit();
  };
  /*!
    \class SM_Phantom_U1
    \brief contains the SM_Phantom_U1 (no surprise)
  */
  /*!
    \var std::string SM_Phantom_U1::m_generator
    The spectrum generator to be used: None.
  */
  /*!
    \fn SM_Phantom_U1::SM_Phantom_U1(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, MODEL::SM. 
    Its parameters (i.e. constants, functions, matrices) are then copied into the newly 
    initialized pointers to the searchable maps. 
  */
  /*!
    \fn void SM_Phantom_U1::FillSpectrum()
    This method is called by the constructor and selects the model framework. 
    There are four parameters (in addition to the SM ones):
    - \f$\tan\beta\f$, the ratio \f$v_2/v_1\f$ of the two vevs of the Higgs doublet
      (\f$v_2\f$) and the Higgs singlet (\f$v_1\f$).  
    - \f$\theta\f$, the mixing angle of the two physical Higgs bosons, \f$\theta\to 0\f$ refers
      to the pure Standard Model.
    - \f$m_{H_1,2}\f$, the mass of the two Higgs bosons.
    - \f$m_{Z'}\f$, the mass of the additional \f$U(1)\f$ \f$Z'\f$ boson.
    - \f$g'_1\f$, the gauge coupling of the additional \f$U(1)\f$.
  */
}

#endif


