#ifndef Input_Output_Handler_H
#define Input_Output_Handler_H

#include "CXXFLAGS_PACKAGES.H"
#include "Message.H"
#ifdef USING__GZIP
#include "Gzip_Stream.H"
#endif

#include <map>
#include <vector>

namespace ATOOLS {
  class Blob_List;
}

namespace SHERPA {

  class HepEvt_Interface;
#ifdef USING__CLHEP
  class HepMC_Interface;
#endif
#ifdef USING__HEPMC2
  class HepMC2_Interface;
#endif

  struct iotype {
    enum code {
      Unknown  = 0,
      Sherpa   = 1,
      HepMC    = 2,
      OldHepMC = 4,
      HepEvt   = 8,
      D0HepEvt = 16,
      HepMC2   = 32,
      size     = 64
    };
  };// end of struct iotype
  
  const iotype::code operator|(const iotype::code code1,const iotype::code code2);  
  const iotype::code operator&(const iotype::code code1,const iotype::code code2);  

  class NameStream {
  public:
    std::string basicfilename, fileextension;
#ifdef USING__GZIP
    ATOOLS::ogzstream outstream;
#else
    std::ofstream     outstream;
#endif
    NameStream(std::string bfname,std::string fext,int precision) :
      basicfilename(bfname), fileextension(fext)
    {
#ifdef USING__GZIP
      fileextension += ".gz";
#endif
      outstream.open((basicfilename+fileextension).c_str());
      outstream.precision(precision);
      if (!outstream.good()) { 
	msg_Error()<<"ERROR in Input_Output_Handler."<<std::endl
		   <<"   Could not open event file "<<(basicfilename+fileextension)<<"."<<std::endl
		   <<"   Will abort the run."<<std::endl;
	abort();
      }
    }
  };

  class Input_Output_Handler {
  private :
    bool            m_on;
    int             m_io,m_precision;
    iotype::code    m_outtype, m_screenout, m_intype;

    ATOOLS::Blob_List        * p_blobs;        
#ifdef USING__CLHEP
    HepMC_Interface          * p_hepmc;
#endif
#ifdef USING__HEPMC2
    HepMC2_Interface         * p_hepmc2;
#endif
    HepEvt_Interface * p_hepevt;
    std::ifstream            * p_instream;
    std::string                m_path, m_filename, m_file;
    long int                   m_filesize, m_evtnumber, m_evtcount;

    std::map<iotype::code,NameStream *> m_outmap;

    bool InitialiseInput(const std::vector<std::string> &);
    bool InitialiseOutput(const std::string, const std::vector<std::string> &);

    void ResetInterfaces();
    void SherpaOutput(std::ostream & outstream,ATOOLS::Blob_List *const blobs,const double=1.);
    bool SherpaInput(ATOOLS::Blob_List *const blobs);
  public :

    // constructors
    Input_Output_Handler(const std::string,
			 const std::vector<std::string> &outfiles,
			 const std::vector<std::string> &infiles,
			 const std::string _path=std::string("./"), 
			 const int _filesize=1000, const int precision=6);

    // destructor
    ~Input_Output_Handler();

    // member functions
    void AddOutputMode(const iotype::code);
    void AddInputMode(const iotype::code);

    bool OutputToFormat(ATOOLS::Blob_List * const,const double=1.);
    bool InputFromFormat(ATOOLS::Blob_List * const);

    void PrintEvent(ATOOLS::Blob_List *const);
    // inline functions
    inline bool OutputOn() const { return (m_io&1 || m_io&2); }

#ifdef USING__CLHEP
    inline HepMC_Interface *const GetHepMCInterface() const { return p_hepmc; }
#endif
#ifdef USING__HEPMC2
    inline HepMC2_Interface *const GetHepMC2Interface() const { return p_hepmc2; }
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
