#ifndef Initialization_Handler_h
#define Initialization_Handler_h

#include "Model_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Matrix_Element_Handler.H"
#include "Hard_Decay_Handler.H"
#include "Shower_Handler.H"
#include "Beam_Remnant_Handler.H"
#include "Fragmentation_Handler.H"
#include "Hadron_Decay_Handler.H"
#include "MI_Handler.H"
#include "Analysis_Handler.H"
#include "Lund_Interface.H"
#include "Event_Reader.H"
#include "Data_Read.H"
#ifdef USING__MCatNLO
#include "Herwig_Interface.H"
#include "MCatNLO_Interface.H"
#endif
#ifdef USING__Ahadic
#include "Ahadic.H"
#endif

#include <map>
#include <string>
#include <iostream> 

namespace SHERPA {

  struct isr {

    enum id {
      none            =  0,
      hard_process    =  1,
      hard_subprocess =  2,
      unknown         = 99
    };

  };// end of struct isr

  typedef std::map<isr::id,PDF::ISR_Handler*> ISR_Handler_Map;

  class IO_Handler;

  class Initialization_Handler {
  private :
    std::string m_path,m_file,m_evtfile;
    std::string m_modeldat, m_beamdat, m_isrdat[2], m_medat, m_midat, m_decaydat, m_showerdat;
    std::string m_beamremnantdat, m_fragmentationdat, m_hadrondecaysdat, m_analysisdat;

    double              m_bunch_splimits[2];
    ATOOLS::Flavour     m_bunch_particles[2];

    ATOOLS::Data_Read * p_dataread;

    int m_mode;
    std::string m_scan_variable;
    double m_scan_begin, m_scan_end, m_scan_value;
    int m_scan_nsteps;
    int m_scan_istep;
    ATOOLS::Parameter_Map m_options;

    MODEL::Model_Base          * p_model;
    BEAM::Beam_Spectra_Handler * p_beamspectra;
    ISR_Handler_Map              m_isrhandlers;
    MEHandlersMap                m_mehandlers;
    Hard_Decay_Handler         * p_harddecays;
    Shower_Handler             * p_showerhandler;
    Beam_Remnant_Handler       * p_beamremnants;
    Fragmentation_Handler      * p_fragmentation;
    HDHandlersMap                m_hdhandlers;
    MI_Handler                 * p_mihandler; 
    IO_Handler                 * p_iohandler;
    Lund_Interface             * p_pythia;
    Event_Reader               * p_evtreader;
#ifdef USING__MCatNLO
    Herwig_Interface           * p_herwig;
    MCatNLO_Interface          * p_mcatnlo;
#endif
#ifdef USING__Ahadic
    AHADIC::Ahadic             * p_ahadic;
#endif
    ANALYSIS::Analysis_Handler * p_analysis;
  private :
    bool InitializeTheIO();
    bool InitializeTheExternalMC();
    bool InitializeTheModel();
    bool InitializeTheBeams();
    bool InitializeThePDFs();
    bool InitializeTheMatrixElements();
    bool InitializeTheUnderlyingEvents();
    bool InitializeTheHardDecays();
    bool InitializeTheShowers();
    bool InitializeTheBeamRemnants();
    bool InitializeTheFragmentation();
    bool InitializeTheHadronDecays();
    bool InitializeTheAnalyses();

    int  ExtractCommandLineParameters(int argc,char * argv[]);
    void CheckFlagConsistency();
    bool CheckBeamISRConsistency();
    void SetParameter(int nr);

  public :

    Initialization_Handler(int argc,char * argv[]);
    Initialization_Handler(std::string=std::string("./"),
			   std::string=std::string("Run.dat"));
    ~Initialization_Handler();

    bool InitializeTheFramework(int =-1);
    bool CalculateTheHardProcesses();

    const std::string Path()    const { return m_path; }
    const std::string File()    const { return m_file; }
    const std::string EvtFile() const { return m_evtfile; }

    Matrix_Element_Handler      * const GetMatrixElementHandler(std::string);
    Hadron_Decay_Handler        * const GetHadronDecayHandler(std::string);
    //ANALYSIS::Analysis_Handler  * const GetAnalysisHandler(std::string);
    inline PDF::ISR_Handler     * const GetISRHandler(const isr::id id) { return m_isrhandlers[id]; }

    MEHandlersMap               * const GetMatrixElementHandlers()       { return &m_mehandlers;     }
    ANALYSIS::Analysis_Handler  * const GetSampleAnalysis()        const { return p_analysis;        }
    MODEL::Model_Base           * const GetModel()                 const { return p_model;           }
    ISR_Handler_Map             * const GetISRHandlers()                 { return &m_isrhandlers;    }
    BEAM::Beam_Spectra_Handler  * const GetBeamSpectraHandler()    const { return p_beamspectra;     }
    Hard_Decay_Handler          * const GetHardDecayHandler()      const { return p_harddecays;      }
    Shower_Handler              * const GetShowerHandler()         const { return p_showerhandler;   }
    Beam_Remnant_Handler        * const GetBeamRemnantHandler()    const { return p_beamremnants;    }
    Fragmentation_Handler       * const GetFragmentationHandler()  const { return p_fragmentation;   }
    HDHandlersMap               * const GetHadronDecayHandlers()         { return &m_hdhandlers;     }
    ATOOLS::Flavour             * const GetBunches()                     { return m_bunch_particles; }
    MI_Handler                  * const GetMIHandler()             const { return p_mihandler;       }
    IO_Handler                  * const GetIOHandler()             const { return p_iohandler;       }
    Lund_Interface              * const GetPythiaInterface()       const { return p_pythia;          }
    Event_Reader                * const GetEventReader()           const { return p_evtreader;       }
#ifdef USING__MCatNLO
    Herwig_Interface            * const GetHerwigInterface()       const { return p_herwig;          }
    MCatNLO_Interface           * const GetMCatNLOInterface()      const { return p_mcatnlo;         }
#endif
#ifdef USING__Ahadic
    AHADIC::Ahadic              * const GetAhadic()                const { return p_ahadic;          }
#endif

    int    const Mode()           const { return m_mode; };
    int    const NumberOfSteps()  const { return m_scan_nsteps; }
    double const ParameterValue() const { return m_scan_value; }

    int    UpdateParameters();
  };

  /*!
    \file 
    \brief Contains the class SHERPA::Initialization_Handler
  */                                  

  /*!
    \class Initialization_Handler
    \brief It is in charge of initialising and deleting all physics modules used during a run.
    
    This class takes care of the proper initialization of all physics modules to be initialized
    during the set up of the event generator. In particular, these contain
    - the model which contains all physics parameters (like coupling constants, 
      mixing matrices, etc.)
    - beam spectra (which is particularly important for physics processes involving for instance
      laser back scattering of photons off electrons)
    - the pdf's through the ISR_Handler,
    - the matrix elements, if needed they are also evaluated to yield total cross sections,
    - the parton shower in the initial and final state,
    - the handling of the beam remnants,
    - the fragmentation which at the moment is dealt with through an interface to the Lund string
      of Pythia and the subsequent hadron decays by the same program.

    These modules are later on taken by the main steering class, Sherpa, and put into 
    corresponding Event_Phase_Handlers. After termination of the run, they are all deleted.
  */
  /*!
    \var std::string Initialization_Handler::m_path
    The path where the steering input data and all other input data files are to be found.
    Also, the particle data are read in through the files "Particle.dat" and "Hadron.dat" 
    in the same path. The default for the path is "./".
  */
  /*!
    \var std::string Initialization_Handler::m_file
    The name of the steering input file, the default is "Run.dat". Here mainly the names of all
    other input files if differing from the default values are to be given.
  */
  /*!
    \var std::string Initialization_Handler::m_modeldat
    The name of the model data file. Apart from the name of the model to be initialized also 
    parameters like, e.g. coupling constants are gathered there. See the corresponding specific
    models for more information about which data are expected. The default is "Model.dat".
  */
  /*!
    \var std::string Initialization_Handler::m_beamdat
    The file with the data specifying the beam treatment. As a minimal requirement, beam particles
    and energies should be given as well as information about which type of beam spectrum
    the user wishes to use. The default is "Beam.dat".
  */
  /*!
    \var std::string Initialization_Handler::m_isrdat
    The file with the data specifying the pdf's. It is important that the particles named in this
    file match with the beam treatment, otherwise the run terminates with a corresponding
    error message. The default is "ISR.dat".
  */
  /*!
    \var std::string Initialization_Handler::m_medat
    The file containing information about which matrix elements are to be used (either some toy
    2->2 processes as internal ones or the full wealth of AMEGIC is available here) and which 
    process(es) is (are) to be calculated. The default for this is "ME.dat".
  */
  /*!
    \var std::string Initialization_Handler::m_showerdat
    Parameters and switches defining the shower handling are to be found here. The default is
    "Shower.dat".
  */
  /*!
    \var std::string Initialization_Handler::m_beamremnantdat
    Information about the beam remnant treatment are to be found in thsi file, defautl for it
    is - again - "Beam.dat".
  */
  /*!
    \var std::string Initialization_Handler::m_fragmentationdat
    The parameters and switches of the fragmentation model to be used are specified in the file
    defined by this word string. The default is "Fragmentation.dat".
  */
  /*!
    \var std::string Initialization_Handler::m_hadrondecaysdat
    For hadron decays which are anticipated, the file with this name will be responsible.
    The default for it is "Fragmentation.dat".   
  */
  /*!
    \var double Initialization_Handler::m_bunch_splimits[2]
    The boundaries for the s' after ISR. These paraermters are read in in the Initialization_Handler
    and passed to the ISR_Handler. They are given in units of the total energy squared under
    consideration.
  */
  /*!
    \var ATOOLS::Flavour Initialization_Handler::m_bunch_particles[2]
    The bunch particles, i.e. the particles that are resolved into their constituents by the pdf's.
  */
  /*!
    \var ATOOLS::Data_Read * Initialization_Handler::p_dataread
    A pointer to the data read object used all over the class. 
  */
  /*!
    \var MODEL::Model_Base * Initialization_Handler::p_model
    A pointer to the model that is initialized by the Initialization_Handler.
  */
  /*!
    \var BEAM::Beam_Spectra_Handler * Initialization_Handler::p_beamspectra
    A pointer to the beam spectra that are initialized by the Initialization_Handler.
  */
  /*!
    \var PDF::ISR_Handler * Initialization_Handler::p_isrhandler
    A pointer to the ISR handling (the pdf's) that is initialized by the Initialization_Handler.
  */
  /*!
    \var Matrix_Element_Handler * Initialization_Handler::p_mehandler
    A pointer to the abstract wrapper around the matrix elements (the pdf's) that are initialized 
    by the Initialization_Handler.
  */
  /*!
    \var Shower_Handler * Initialization_Handler::p_showerhandler
    A pointer to the abstract wrapper around the shower module. It is initialized by the 
    Initialization_Handler. 
  */
  /*!
    \var Beam_Remnant_Handler * Initialization_Handler::p_beamremnants
    A pointer to the beam remnant handling used thoughout the run. 
  */
  /*!
    \var Fragmentation_Handler * Initialization_Handler::p_fragmentation
    A pointer to the wrapper around the fragmentation code that will be used during the run.
  */
  /*!
    \var Hadron_Decay_Handler * Initialization_Handler::p_hadrondecays
    A pointer to the wrapper around the hadron decay code that will be used during the run.
  */
  /*!
    \fn Initialization_Handler::Initialization_Handler(std::string=std::string("./"),
                                                       std::string=std::string("Run.dat"))
    In the constructor, the data read object is instantiated and all file names to be used
    for the set up are read in.
  */
  /*!
    \fn Initialization_Handler::~Initialization_Handler()
    The destructor deletes all physdics modules that have been used during the run.
  */
  /*!
    \fn bool Initialization_Handler::InitializeTheFramework()
    This method first reads in all particles, before subsequently the model, the beams and the pdf's
    are set up. After this has been done, consistency of the bunch particles and the beam 
    handling is checked for. If successul, the matrix elements, showers, beam remanants, 
    fragmentation and hadron decays are instantiated. If everything worked out, true is returned
    and some first corss section evaluation will start.
  */
  /*!
    \fn bool Initialization_Handler::CalculateTheHardProcesses()
    Returns true if the calculation of total cross sections by the matrix elements
    of choice was successful, otherwise, flase will be returned.
  */
}


#endif
