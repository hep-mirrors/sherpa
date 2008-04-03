#ifndef MUED_h
#define MUED_h

#include "Model_Base.H"

namespace MODEL {
  class MUED : public Model_Base {
  private :
    void FillSpectrum();
    void RunSpectrumGenerator();
  public :
    MUED(std::string,std::string);
    ~MUED();
    bool FillDecay(ATOOLS::Decay_Table *) { return true; }
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
    \fn void MUED::FillSpectrum()
    This method is called by the constructor and set all the model parameters, to do
    so it subsequently calls MUED::RunSpectrumGenerator().
  */
  /*!
    \fn bool MUED::RunSpectrumGenerator()
    In this method the actual spectrum is obtained employing the class MUED_Spectrum.
  */
}

#endif


