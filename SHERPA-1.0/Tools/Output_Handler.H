#ifndef Output_Handler_H
#define Output_Handler_H

#ifdef _USE_HEPMC_
#include "HepMC_Interface.H"
#endif

#include "HepEvt_Interface.H"
#include "Blob_List.H"

namespace SHERPA {

    struct iotype {

      enum code {
	Unknown = 0,
	Sherpa  = 1,
	HepMC   = 2,
	HepEvt  = 4,
	size    = 8
      };

    };// end of struct iotype

    const iotype::code operator|(const iotype::code code1,const iotype::code code2);  
    const iotype::code operator&(const iotype::code code1,const iotype::code code2);  

  class Output_Handler {
  private :
    bool            m_on;
    int             m_io;
    iotype::code    m_outtype, m_intype;

#ifdef _USE_HEPMC_
    HepMC_Interface          * p_hepmc;
    HepMC::GenEvent          * p_event;
#endif
    ATOOLS::HepEvt_Interface * p_hepevt;
    std::ofstream              m_outstream;
    std::ifstream            * p_instream;
    std::string                m_path, m_filename, m_file;
    long int                   m_filesize, m_evtnumber, m_evtcount;

    void SherpaOutput(ATOOLS::Blob_List *const blobs,const double=1.);
    bool SherpaInput(ATOOLS::Blob_List *const blobs);
  public :

    // constructors
    Output_Handler();
    Output_Handler(const std::vector<std::string> &outfiles,
		   const std::vector<std::string> &infiles,
		   const std::string=std::string("."));

    // destructor
    ~Output_Handler();

    // member functions
    void AddOutputMode(const iotype::code);
    void AddInputMode(const iotype::code);

    bool OutputToFormat(ATOOLS::Blob_List * const,const double=1.);
    bool InputFromFormat(ATOOLS::Blob_List * const);
    
    // inline functions
    inline bool OutputOn() const { return m_io&1; }
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
