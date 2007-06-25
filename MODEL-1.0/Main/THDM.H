#ifndef THDM_h
#define THDM_h

#include "Model_Base.H"

namespace MODEL {
  class THDM : public Model_Base {
  private :
    void ReadInFile();
  public :
    THDM(std::string,std::string);
    ~THDM()                     {}
    bool RunSpectrumGenerator() {}
    void FillMasses();
    void FillWidths() {}
    bool FillDecay(ATOOLS::Decay_Table *) { return 0; }
  };
  /*!
    \class THDM
    \brief contains the THDM (no surprise)
  */
  /*!
    \var std::string THDM::m_generator
    The spectrum generator to be used: Internal.
  */
  /*!
    \fn THDM::THDM(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, MODEL::SM. 
    Its parameters (i.e. constants, functions, matrices) are then copied into the newly 
    initialized pointers to the searchable maps. 
  */
  /*!
    \fn bool THDM::RunSpectrumGenerator()
    In this method the spectrum generator is specified and initialised to run.
  */
  /*!
    \fn void THDM::ReadInFile()
    This method is called by the constructor and selects the model framework.
  */
  /*!
    \fn void THDM::FillMasses()
    \todo This method should be implemented.
  */
  /*!
    \fn void THDM::FillWidths()
    This method should be implemented. A good way would be to calculate couplings in vertices
    from decay widths to guarantee internal consistency.
  */
}

#endif


