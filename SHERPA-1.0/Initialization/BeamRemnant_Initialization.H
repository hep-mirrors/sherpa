#ifndef BeamRemnant_Initialization_h
#define BeamRemnant_Initialization_h

#include "Model_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Beam_Remnant_Handler.H"

#include <string>
#include <iostream> 

namespace SHERPA {
  class BeamRemnant_Initialization {
  private :
    std::string                  m_dir,m_file;
    AORGTOOLS::Data_Read       * p_dataread;

    MODEL::Model_Base          * p_model;
    ISR::ISR_Handler           * p_isr;
    BEAM::Beam_Spectra_Handler * p_beam;
    Beam_Remnant_Handler       * p_beamremants;
  private :
    void ReadInFile();
    bool InitializeTheBeamRemnantHandler();
  public :
    BeamRemnant_Initialization(std::string,std::string,MODEL::Model_Base *,
			       ISR::ISR_Handler *,BEAM::Beam_Spectra_Handler *);
    ~BeamRemnant_Initialization();

    void Output();
  };
}


#endif
