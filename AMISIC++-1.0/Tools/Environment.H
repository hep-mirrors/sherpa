#ifndef Environmen_H
#define Environmen_H

#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Model_Base.H"

namespace AMISIC {

  class Environment {
  private:

    std::string m_path, m_file;
    std::string m_modeldat, m_beamdat, m_isrdat;

    BEAM::Beam_Spectra_Handler *p_beamhandler;
    PDF::ISR_Handler           *p_isrhandler;
    MODEL::Model_Base          *p_model;

    ATOOLS::Flavour m_bunch[2], m_beam[2];

  private :

    bool InitializeTheModel();
    bool InitializeTheBeams();
    bool InitializeThePDFs();

  public :

    // constructor
    Environment(const std::string,const std::string);

    // destructor
    ~Environment();

    // member functions
    bool InitializeTheEnvironment();

    MODEL::Model_Base          *Model()       { return p_model;       }
    BEAM::Beam_Spectra_Handler *BeamHandler() { return p_beamhandler; }
    PDF::ISR_Handler           *ISRHandler()  { return p_isrhandler;  }

  };// end of class Environment

  //////////////////
  // Doxygen part //
  //////////////////

  /*!
    \file Environment.H
    \brief Declares the class Environment
  */

  /*!
    \class Environment
    \brief Allows running Amisic as a stand-alone 
    multiple interaction generator
  */

}// end of namespace AMISIC

#endif
