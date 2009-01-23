#ifndef Input_Output_Handler_H
#define Input_Output_Handler_H

#include "CXXFLAGS_PACKAGES.H"
#include "Message.H"

#include <map>
#include <vector>

namespace ATOOLS {
  class Blob_List;
  class Data_Reader;
}

namespace SHERPA {

#ifdef USING__HEPMC2
  class HepMC2_Interface;
#endif
  class HepEvt_Interface;
  class Output_Base;

  class Input_Output_Handler {
  private :
    std::map<std::string, Output_Base *> m_outmap;
    int m_filesize, m_evtnumber, m_evtcount;
#ifdef USING__HEPMC2
    HepMC2_Interface         * p_hepmc2;
#endif
    HepEvt_Interface         * p_hepevt;
    std::string m_outmode;

    bool InitialiseOutput(ATOOLS::Data_Reader* dr);
    void ResetInterfaces();
  public :

    // constructors
    Input_Output_Handler(ATOOLS::Data_Reader* dr);

    // destructor
    ~Input_Output_Handler();

    bool OutputToFormat(ATOOLS::Blob_List * const);
    bool InputFromFormat(ATOOLS::Blob_List * const) { return false; }

    void PrintEvent(ATOOLS::Blob_List *const);

#ifdef USING__HEPMC2
    inline HepMC2_Interface *const GetHepMC2Interface() const { 
      return p_hepmc2; 
    }
#endif
  };

  /*!
    \class Input_Output_Handler
    \brief Handles the output of single events to an output format.

    This class steers the output of single events to a specific output selected
    from outside. So far, two formats are available (and mildly tested):
    HepMC through the HepMC_Interface and HepEvt through the HepEvt_Interface.
  */
  /*!
    \var HepMC_Interface * Input_Output_Handler::p_hepmc
  */
  /*!
    \var HepMC::GenEvent * Input_Output_Handler::p_event
  */
  /*!
    \var bool Input_Output_Handler::m_active
  */
  /*!
    \var int Input_Output_Handler::m_type
  */
  /*!
    \fn Input_Output_Handler::Input_Output_Handler(int type=1)
  */
  /*!
    \fn Input_Output_Handler::~Input_Output_Handler()
  */
  /*!
    \fn void Input_Output_Handler::OutputToFormat(ATOOLS::Blob_List *)
  */
  /*!
    \fn bool Input_Output_Handler::Active()
  */
}
#endif
