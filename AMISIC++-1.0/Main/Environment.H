#ifndef Environmen_h
#define Environmen_h

#include "Data_Read.H"
#include "Run_Parameter.H"

#include "Model_Base.H"
#include "Beam_Spectra_Handler.H"
#include "Beam_Base.H"
#include "ISR_Handler.H"
#include "ISR_Base.H"

#include "Flavour.H"

#include <string>
#include <iostream> 

namespace AMEGIC {
  class Environment {
  private :
    std::string m_path,m_file;
    std::string m_modeldat,m_beamdat, m_isrdat, m_medat;

    ATOOLS::Data_Read       * p_dataread;

    MODEL::Model_Base          * p_model;

    BEAM::Beam_Spectra_Handler * p_beamspectra;
    ATOOLS::Flavour           m_beam_particles[2];

    PDF::ISR_Handler           * p_isrhandler;
    ATOOLS::Flavour           m_bunch_particles[2];
    double                       m_bunch_splimits[2];

  private :
    bool InitializeTheModel();
    bool InitializeTheBeams();
    bool InitializeThePDFs();
    //bool InitializeTheMatrixElements();
  public :
    Environment(std::string,std::string);
    void InitializeTheEnvironment();
    ~Environment();

    MODEL::Model_Base          * Model()              { return p_model; }
    BEAM::Beam_Spectra_Handler * BeamSpectraHandler() { return p_beamspectra; }
    PDF::ISR_Handler           * ISRHandler()         { return p_isrhandler; }

    std::string                  MEFile()             { return m_medat; }
  };
}


#endif
