#ifndef Fourth_Generation_Leptons_h
#define Fourth_Generation_Leptons_h

#include "Model_Base.H"

namespace MODEL {
  class Fourth_Generation_Leptons : public Model_Base {
  private :
    Model_Base *p_sm;
    void ReadInFile();
  public :
    Fourth_Generation_Leptons(std::string,std::string);
    ~Fourth_Generation_Leptons();
    bool RunSpectrumGenerator() { return true; }
    void FillMasses();
    void FillWidths() {}
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
    \fn void Fourth_Generation_Leptons::ReadInFile()
    This method is called by the constructor and selects the model framework.
  */
  /*!
    \fn void Fourth_Generation_Leptons::FillMasses()
    \todo This method should be implemented.
  */
  /*!
    \fn void Fourth_Generation_Leptons::FillWidths()
    This method should be implemented. A good way would be to calculate couplings in vertices
    from decay widths to guarantee internal consistency.
  */
}

#endif


