#ifndef MODEL_Main_MSSM_H
#define MODEL_Main_MSSM_H

#include "MODEL/Main/Model_Base.H"

namespace MODEL {
  class MSSM : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
    void RunSpectrumGenerator();
  public :
    MSSM(std::string,std::string,bool);
    ~MSSM();
  };
  /*!
    \class MSSM
    \brief contains the MSSM (no surprise)
  */
  /*!
    \fn MSSM::MSSM(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, MODEL::SM. 
    Its parameters (i.e. constants, functions, matrices) are then copied into the newly 
    initialized pointers to the searchable maps. 
  */
  /*!
    \fn bool MSSM::RunSpectrumGenerator()
    This method creates and calls the LesHouches_Interface class to read in the SUSY
    spectrum information.
  */
  /*!
    \fn void MSSM::FillSpectrum()
    This method is responsable for setting all the SUSY parameters. For this 
    SUSY Les Houches Accord conform input files are read through the LesHouches_Interface
    class that is build up through the method MSSM::RunSpectrumGenerator()
  */
}

#endif


