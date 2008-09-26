#ifndef SM_AxiGluon_h
#define SM_AxiGluon_h

#include "Model_Base.H"
#include "MyComplex.H"

namespace MODEL {
  class SM_AxiGluon : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
    bool SanityChecks();
  public :
    SM_AxiGluon(std::string,std::string,bool);
    ~SM_AxiGluon();
    bool FillDecay(ATOOLS::Decay_Table *) { return 0; }
  };
  /*!
    \class SM_AxiGluon
    \brief contains the SM_AxiGluon (no surprise)
  */
  /*!
    \var std::string SM_AxiGluon::m_generator
    The spectrum generator to be used: None.
  */
  /*!
    \fn SM_AxiGluon::SM_AxiGluon(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, MODEL::SM. 
    Its parameters (i.e. constants, functions, matrices) are then copied into the newly 
    initialized pointers to the searchable maps. 
  */
  /*!
    \fn void SM_AxiGluon::FillSpectrum()
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


