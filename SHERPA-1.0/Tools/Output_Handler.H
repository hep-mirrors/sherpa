#ifndef Output_Handler_h
#define Output_Handler_h

#ifdef _USE_HEPMC_
#include "HepMC_Interface.H"
#endif

#include "HepEvt_Interface.H"
#include "Blob_List.H"

namespace SHERPA {
  class Output_Handler {
  private :
#ifdef _USE_HEPMC_
    HepMC_Interface  * p_hepmc;
    HepMC::GenEvent  * p_event;
#endif
    HepEvt_Interface * p_hepevt;
    bool               m_active;
    int                m_type;
  public :
    Output_Handler(int type=1);
    ~Output_Handler();
    void OutputToFormat(ATOOLS::Blob_List *);
    bool Active() { return m_active; }
  };

  /*!
    \class Output_Handler
    \brief Handles the output of single events to an output format.

    This class steers the output of single events to a specific output selected
    from outside. So far, two formats are available (and mildly tested):
    HepMC through the HepMC_Interface and HepEvt through the HepEvt_Interface.
  */
  /*!
    \var HepMC_Interface * Output_Handler::p_hepmc
  */
  /*!
    \var HepMC::GenEvent * Output_Handler::p_event
  */
  /*!
    \var bool Output_Handler::m_active
  */
  /*!
    \var int Output_Handler::m_type
  */
  /*!
    \fn Output_Handler::Output_Handler(int type=1)
  */
  /*!
    \fn Output_Handler::~Output_Handler()
  */
  /*!
    \fn void Output_Handler::OutputToFormat(ATOOLS::Blob_List *)
  */
  /*!
    \fn bool Output_Handler::Active()
  */
}
#endif
