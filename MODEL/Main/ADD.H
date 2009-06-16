#ifndef MODEL_Main_ADD_H
#define MODEL_Main_ADD_H

#include "MODEL/Main/Model_Base.H"

namespace MODEL {
  class ADD : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
  public :
    ADD(std::string,std::string,bool);
    ~ADD();
  };
  /*!
    \class ADD
    \brief A model for extra dimensions accroding to ADD.
  */
  /*!
    \fn ADD::ADD(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, MODEL::SM. 
    Its parameters (i.e. constants, functions, matrices) are then copied into the newly 
    initialized pointers to the searchable maps. 
  */
  /*!
    \fn bool ADD::RunSpectrumGenerator()
    In this model, no spectrum generator is needed so far.
  */
  /*!
    \fn void ADD::FillSpectrum()
    This method is called by the constructor and defines the model parameters.
    The following default values apply:
    -f[G_N = 6.707\cdot 10^{-39},\; 
       N_{\rm E.D.} = 2,\;
       M_S = r = 0\;.\f]
    From them, all other parameters are directly calculated inside this method.
  */
}

#endif


