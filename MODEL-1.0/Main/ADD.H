#ifndef ADD_h
#define ADD_h

#include "Model_Base.H"

namespace MODEL {
  class ADD : public Model_Base {
  private :
    void   ReadInFile();
  public :
    ADD(std::string,std::string);
    ~ADD()                                {}
    bool RunSpectrumGenerator()           { return 1; }
    void FillMasses()                     {}
    void FillWidths()                     {}
    bool FillDecay(ATOOLS::Decay_Table *) { return 0; }
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
    \fn void ADD::ReadInFile()
    This method is called by the constructor and defiens the model parameters.
    The following default values apply:
    -f[G_N = 6.707\cdot 10^{-39},\; 
       N_{\rm E.D.} = 2,\;
       M_S = r = 0\;.\f]
    From them, all other parameters are directly calculated inside this method.
  */
}

#endif


