#ifndef SHERPA_Initialization_Initialization_Handler_H
#define SHERPA_Initialization_Initialization_Handler_H

#include "ATOOLS/Org/CXXFLAGS.H"

#include "MODEL/Main/Model_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"
#include "SHERPA/SoftPhysics/Fragmentation_Handler.H"
#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"
#include "SHERPA/LundTools/Lund_Interface.H"
#include "SHERPA/Tools/Event_Reader.H"
#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/Data_Reader.H"

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

  class Input_Output_Handler;

  class Initialization_Handler: public ATOOLS::Terminator_Object {
  private :
    std::string m_path,m_file,m_evtfile,m_pdflib;
    std::string m_modeldat, m_beamdat, m_isrdat[2], m_medat, m_midat, m_showerdat;
    std::string m_beamremnantdat, m_fragmentationdat, m_hadrondecaysdat, m_analysisdat;
    std::string m_softphotonsdat;
    std::string m_processesdat;
    std::string m_selectordat;

    double              m_bunch_splimits[2];
    ATOOLS::Flavour     m_bunch_particles[2];

    ATOOLS::Data_Reader *p_dataread;

    int m_mode;
    bool m_spincorrelations, m_savestatus;

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
    Soft_Photon_Handler        * p_softphotons;
    Input_Output_Handler       * p_iohandler;
    Lund_Interface             * p_pythia;
    Event_Reader               * p_evtreader;
    Analysis_Map                 m_analyses;
    void ShowParameterSyntax();
    void LoadLibraries() const;
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
    bool InitializeTheSoftPhotons();
    bool InitializeTheAnalyses();

    bool ExtractValArg
    (std::vector<std::string> &args,
     std::vector<std::string>::iterator &it,
     const std::string &arg,const std::string &tag,
     const std::string &def=std::string("")) const;
    int  ExtractCommandLineParameters(int argc,char * argv[]);
    bool CheckBeamISRConsistency();
    void SetScaleFactors();
    void PrepareTerminate();
    void SetFileNames();

  public :

    Initialization_Handler(int argc,char * argv[]);
    ~Initialization_Handler();

    bool InitializeTheFramework(int =-1);
    bool CalculateTheHardProcesses();

    const std::string Path()    const { return m_path; }
    const std::string File()    const { return m_file; }
    const std::string EvtFile() const { return m_evtfile; }

    Matrix_Element_Handler      * const GetMatrixElementHandler(std::string);
    Hadron_Decay_Handler        * const GetHadronDecayHandler(std::string);
    inline PDF::ISR_Handler     * const GetISRHandler(const isr::id id) { return m_isrhandlers[id]; }
    MEHandlersMap               * const GetMatrixElementHandlers()       { return &m_mehandlers;     }
    Analysis_Map                * const GetAnalyses()                    { return &m_analyses;        }
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
    Soft_Photon_Handler         * const GetSoftPhotonHandler()     const { return p_softphotons;     }
    Input_Output_Handler        * const GetIOHandler()             const { return p_iohandler;       }
    Lund_Interface              * const GetPythiaInterface()       const { return p_pythia;          }
    Event_Reader                * const GetEventReader()           const { return p_evtreader;       }

    int    const Mode()           const { return m_mode; };

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
    \var HDHandlersMap Initialization_Handler::m_hdhandlers
    This is map with pointers to the corresponding hadron decay handlers
    <TABLE>
    <TR> <TD>1st column</TD> <TD>2nd column</TD> </TR>
    <TR> <TD>Lund</TD> <TD>pointer on Pythia</TD> </TR>
    <TR> <TD>Hadrons</TD> <TD>pointer on HADRONS</TD> </TR>
    </TABLE>
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

  /*!
    \fn Initialization_Handler::InitializeTheHadronDecays()
    This method reads \c Fragmentation.dat, especially the desired decay model.
     
    If the model is \c Hadrons it calls the HADRONS constructor for initialisation
    and the Hadron_Decay_Handler::EraseTreated method
    to remove the particles that HADRONS can cope with from the unstable hadron list. If then
    there are still some particles in the list, i.e. those that HADRONS can't deal with, it needs to
    initialise Pythia as well.

    If, on the other hand, the model is \c Lund, it calls the Pythia for initialisation.
    */
  

}


#endif
