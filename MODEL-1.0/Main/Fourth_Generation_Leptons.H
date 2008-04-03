#ifndef Fourth_Generation_Leptons_h
#define Fourth_Generation_Leptons_h

#include "Model_Base.H"

namespace MODEL {
  class Fourth_Generation_Leptons : public Model_Base {
  private :
    void FillSpectrum();
  public :
    Fourth_Generation_Leptons(std::string,std::string);
    ~Fourth_Generation_Leptons();
    bool FillDecay(ATOOLS::Decay_Table *) { return 0; }
  };
  /*!
    \class Fourth_Generation_Leptons
    \brief contains the Fourth_Generation_Leptons (no surprise)
  */
  /*!
    \var std::string Fourth_Generation_Leptons::m_generator
    The spectrum generator to be used: Internal.
  */
  /*!
    \fn Fourth_Generation_Leptons::Fourth_Generation_Leptons(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, MODEL::SM. 
    Its parameters (i.e. constants, functions, matrices) are then copied into the newly 
    initialized pointers to the searchable maps. 
  */
  /*!
    \fn bool Fourth_Generation_Leptons::RunSpectrumGenerator()
    In this method the spectrum generator is specified and initialised to run.
  */
  /*!
    \fn void Fourth_Generation_Leptons::FillSpectrum()
    This method is called by the constructor and selects the model framework.
  */
}

#endif


