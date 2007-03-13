#ifndef MUED_h
#define MUED_h

#include "Model_Base.H"

namespace MODEL {
  class MUED : public Model_Base {
  private :
    std::string m_generator;
    void ReadInFile();
  public :
    MUED(std::string,std::string);
    ~MUED()                   {}
    bool RunSpectrumGenerator();
    void FillMasses()                     {}
    void FillWidths()                     {}
    bool FillDecay(ATOOLS::Decay_Table *) { return true; }
    void AddToDataCollector();
  };
  /*!
    \class MUED
    \brief contains the MUED (no surprise)    
  */
  /*!
    \var std::string MUED::m_generator
    The spectrum generator to be used. At the moment, this is all done internally.
  */
  /*!
    \var double M_R
    The size of the extra dimensions.
  */
  /*!
    \fn MUED::MUED(std::string,std::string)
    The constructor initializes the parameters of this class and the pointers to
    the searchable maps.
  */
  /*!
    \fn bool MUED::RunSpectrumGenerator()
    In this method the spectrum generator is specified and initialised to run.
  */
  /*!
    \fn void MUED::ReadInFile()
    This method is called by the constructor.
  */
  /*!
    \fn void MUED::FillMasses()
    \todo This method should be implemented.
  */
  /*!
    \fn void MUED::FillWidths()
    This method should be implemented. A good way would be to calculate couplings in vertices
    from decay widths to guarantee internal consistency.
  */
}

#endif


