#ifndef IO_Handler_H
#define IO_Handler_H

#ifdef USING__CLHEP
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

  class IO_Handler {
  private :
    bool            m_on;
    int             m_io;
    iotype::code    m_outtype, m_intype;

    ATOOLS::Blob_List        * p_blobs;        
#ifdef USING__CLHEP
    HepMC_Interface          * p_hepmc;
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
    IO_Handler(const std::string=std::string("Sherpa"));
    IO_Handler(const std::vector<std::string> &outfiles,
	       const std::vector<std::string> &infiles,
	       const std::string _path=std::string("./"), 
	       const int _filesize=1000);

    // destructor
    ~IO_Handler();

    // member functions
    void AddOutputMode(const iotype::code);
    void AddInputMode(const iotype::code);

    bool OutputToFormat(ATOOLS::Blob_List * const,const double=1.);
    bool InputFromFormat(ATOOLS::Blob_List * const);

    void PrintEvent();
    // inline functions
    inline bool OutputOn() const { return (m_io&1 || m_io&2); }

#ifdef USING__CLHEP
    inline HepMC_Interface *const GetHepMCInterface() const 
    { return p_hepmc; }
#endif

  };

  /*!
    \class IO_Handler
    \brief Handles the output of single events to an output format.

    This class steers the output of single events to a specific output selected
    from outside. So far, two formats are available (and mildly tested):
    HepMC through the HepMC_Interface and HepEvt through the HepEvt_Interface.
  */
  /*!
    \var HepMC_Interface * IO_Handler::p_hepmc
  */
  /*!
    \var HepMC::GenEvent * IO_Handler::p_event
  */
  /*!
    \var bool IO_Handler::m_active
  */
  /*!
    \var int IO_Handler::m_type
  */
  /*!
    \fn IO_Handler::IO_Handler(int type=1)
  */
  /*!
    \fn IO_Handler::~IO_Handler()
  */
  /*!
    \fn void IO_Handler::OutputToFormat(ATOOLS::Blob_List *)
  */
  /*!
    \fn bool IO_Handler::Active()
  */
}
#endif
