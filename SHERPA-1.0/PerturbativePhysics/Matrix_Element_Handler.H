#ifndef Matrix_Element_Handler_H
#define Matrix_Element_Handler_H

#include <map> 
#include <iostream> 

#include "Model_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Particle_List.H"
#include "Particle.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Data_Reader.H"
// #include "Spin_Correlation_Tensor.H"
#include "Amplitude_Tensor.H"

namespace AMEGIC { class Amegic; 
                   class Point;  }

namespace EXTRAXS { class Simple_XS; 
                    class XS_Base;     }

namespace SHERPA {
  typedef std::map<ATOOLS::Flavour,ATOOLS::Flavour> Flavour_Map;

  class Particle_Map {
    bool m_flip_anti;
    bool m_change_flavs;
    Flavour_Map m_flmap;
  public:
    Particle_Map();
    bool Unity();
    void SetFlipAnti();
    void Add(const ATOOLS::Flavour & ,const ATOOLS::Flavour &);
    bool Apply(int n, ATOOLS::Flavour *, ATOOLS::Vec4D *);
  };


  class Matrix_Element_Handler {
  private :

    std::string          m_dir,m_file;
    AMEGIC::Amegic     * p_amegic;
    EXTRAXS::Simple_XS * p_simplexs;
    PDF::ISR_Handler   * p_isr;
    MODEL::Model_Base  * p_model;
    int                  m_mode;
    double               m_weight, m_procweight, m_xsecweight, m_xf1, m_xf2;
    unsigned long          m_ntrial, m_xsecntrial, m_sntrial;

    std::string            m_name,m_signalgenerator,m_readin;
 
    int  m_eventmode;
    int  m_sudakovon, m_apply_hhmf;
    bool m_ini_swaped;
    int  m_nmoms;
    bool m_spincorrelations;

    ATOOLS::Data_Reader    * p_reader;

    ATOOLS::Flavour *  p_flavs;
    ATOOLS::Vec4D   *  p_moms;

    std::vector<Particle_Map> m_particle_maps;

    bool  PrintTotalXSec() ;
    int   InitializeSimpleXS(MODEL::Model_Base *,
			      BEAM::Beam_Spectra_Handler *,PDF::ISR_Handler *);
    int   InitializeAmegic(MODEL::Model_Base *,
			    BEAM::Beam_Spectra_Handler *,PDF::ISR_Handler *);

    void  SetupHHMF();
    bool  ApplyHHMF();
    bool  GetMomentaNFlavours();

    ATOOLS::Blob_Data_Base * UnweightedEvent();
    ATOOLS::Blob_Data_Base * UnweightedSameEvent();

    ATOOLS::Blob_Data_Base * WeightedEvent();
    ATOOLS::Blob_Data_Base * WeightedSameEvent();

    
  public :

    Matrix_Element_Handler();
    Matrix_Element_Handler(std::string,std::string,MODEL::Model_Base *,
			   Matrix_Element_Handler *);
    Matrix_Element_Handler(std::string,std::string,MODEL::Model_Base *,
			   BEAM::Beam_Spectra_Handler *,PDF::ISR_Handler *,
			   Matrix_Element_Handler *);

    ~Matrix_Element_Handler();

    bool InitializeProcesses(BEAM::Beam_Spectra_Handler *,PDF::ISR_Handler *);
    bool InitializeDecays(bool);
    bool AddToDecays(const ATOOLS::Flavour &);
    bool AddToDecays(ATOOLS::Decay_Channel *);
    bool InitializeDecayTables(); 
    bool CalculateWidths();
    bool FillDecayTable(ATOOLS::Decay_Table *,bool); 
    bool CalculateTotalXSecs(int);
    bool PrepareXSecTables();
    bool LookUpXSec(double,bool,std::string);

    bool GenerateOneEvent();
    bool GenerateOneEvent(ATOOLS::Decay_Channel *,double);
    bool GenerateSameEvent();

    double ExpectedEvents() const;
    void   AddEvent(const double xs,const double validxs,const int ncounts);

    void SetAmegic(AMEGIC::Amegic *);
    void SetXS(EXTRAXS::Simple_XS *);

    std::string          SignalGenerator();
    std::string          ProcessName();
    std::string          Name();
    unsigned int         NIn();
    unsigned int         NOut();
    unsigned int         NDecOut();
    unsigned int         MaxJets();
    unsigned int         MinQCDJets();
    unsigned int         MaxQCDJets();
    const ATOOLS::Vec4D      * Momenta();
    const ATOOLS::Vec4D      * DecMomenta();
    const ATOOLS::Flavour    * Flavours();
    int                  InSwaped();

    int                  OrderStrong();
    int                  OrderEWeak();
    unsigned int         NumberOfDiagrams();
    AMEGIC::Amegic     * GetAmegic();
    AMEGIC::Point      * GetDiagram(int);
    EXTRAXS::XS_Base   * GetXS(const int mode=0);
//     ATOOLS::Spin_Correlation_Tensor* GetSpinCorrelations();
    void                 FillAmplitudes(ATOOLS::Amplitude_Tensor*);
    bool                 SpinCorrelations();
    void                 SetSpinCorrelations(bool sc);
    double               Weight();
    double               ProcessWeight();
    double               XSecWeight();
    unsigned long        NumberOfTrials();
    unsigned long        NumberOfXSecTrials();

    PDF::ISR_Handler   * GetISR_Handler()   { return p_isr; }
    MODEL::Model_Base  * GetModel()         { return p_model; }
    int                  UseSudakovWeight() { return m_sudakovon; }
    double               FactorisationScale();

    inline void SetUseSudakovWeight(const int on) { m_sudakovon=on; }
    inline int  EventGenerationMode() const { return m_eventmode; }
    inline void SaveNumberOfTrials() { m_sntrial+=m_ntrial; }
    inline void ResetNumberOfTrials() { m_sntrial=0; }

    inline double XF1() const { return m_xf1; }
    inline double XF2() const { return m_xf2; }

    void SetModel(MODEL::Model_Base *model) { p_model=model; }

  };

  typedef std::map<std::string,Matrix_Element_Handler *> MEHandlersMap;
  typedef MEHandlersMap::iterator                        MEHandlerIter;

  /*!
    \file 
    \brief Contains the class SHERPA::Matrix_Element_Handler
  */
  /*!
    \class Matrix_Element_Handler
    \brief A wrapper around matrix element generators.

    This class provides the generic interface of variuos matrix element generators
    with the rest of the SHERPA framework. As it stands, SHERPA employs their
    matrix elements in the Event_Phase_Handler Signal_Processes and in the
    Event_Phase_Handler Jet_Evolution for the interface with the corresponding
    Shower_Handler. Therefore this handler class provides methods:
    - to initialize the specific matrix element generators selected by the user
      including the initialization of corresponding processes,
    - to calculate their total cross sections within some cuts specified by 
      appropriate selector data (files)
    - to generate weighted and unweighted single events on the parton level
    - to pass more information to other places in the program, like for instance
      into the interface with the shower. This information includes not only 
      the flavours and momenta of the partons but also the Feynman graphs. 
      It is anticipated that in the future methods to extract the colour flow
      will also become available.

    In case the user wants to niterface to another matrix element generator, this 
    handler has to be modified accordingly. Also, an interface to the shower has to 
    be provided to take full advantage of the event generator structure of SHERPA.
  */
  /*!
    \var std::string Matrix_Element_Handler::m_dir
    The directory for the file with the matrix element data.
  */
  /*!
    \var std::string Matrix_Element_Handler::m_file
    The name of the data file with the matrix element specificiations.
  */
  /*!
    \var AMEGIC::Amegic * Matrix_Element_Handler::p_amegic
    A pointer to the AMEGIC object.
  */
  /*!
    \var EXTRAXS::SimpleXSecs * Matrix_Element_Handler::p_simplexs
    A pointer to the list of internal simple \f$2\to 2\f$ processes.
  */
  /*!
    \var int Matrix_Element_Handler::m_mode
    The mode in which the matrix elements are used:
    0 - the initialization did not work out. 
    1 - Amegic is being used.
    2 - The simple internal matrix elements are used.
  */
  /*!
    \var double Matrix_Element_Handler::m_weight
    The weight if the run is in the "weighted events"-mode. The weight is given in
    picobarns.
  */
  /*!
    \var std::string Matrix_Element_Handler::m_name
    The name of the matrix elements used. At the moment, two different matrix element
    objects are available, namely Amegic and some internal simple ones. The respective
    names are "Amegic" and "Simple X-section".
  */
  /*!
    \var std::string Matrix_Element_Handler::m_signalgenerator
    The user's choice of signal generator to be read in through the matrix element data.
    The default is "Amegic", at the moment, also "Internal" for the very limited set of simple 
    processes can be chosen as well.
  */
  /*!
    \var std::string Matrix_Element_Handler::m_readin
    The directory with stored results. This is an extremely dangerous option to use, because
    no real check is done whether the results stored there (and therefore used throughout the run)
    coincide with the actual set-up of the program.
  */
  /*!
    \var bool Matrix_Element_Handler::m_eventmode
    The mode in which events are generated, either "Weighted" or "Unweighted", the latter is the
    default.
  */
  /*!
    \var ATOOLS::Data_Read * Matrix_Element_Handler::p_dataread
    A pointer to the Data_Read object.
  */
  /*!
    \fn bool Matrix_Element_Handler::RescaleJetrates() 
    A method to be filled later. It is meant to rescale jetrates etc..
  */
  /*!
    \fn int Matrix_Element_Handler::InitializeSimpleXS(MODEL::Model_Base *,
                                                       BEAM::Beam_Spectra_Handler *,
                                                       PDF::ISR_Handler *)
    Initializes the user-defined set of simple matrix elements available in  EXTRAXS. More work
    has to be spent there. In principle, the idea is similar to the one for Amegic. First
    the wrapper around hte individual processes is instantiated, then the processes themselves
    are initialized through a call to InitializeProcesses. If this was successful, a 2 is returned, 
    otherwise a 0 is returned, resulting in a controlled end of the run.
  */
  /*!
    \fn int Matrix_Element_Handler::InitializeAmegic(MODEL::Model_Base *,
			                             BEAM::Beam_Spectra_Handler *,
                                                     PDF::ISR_Handler *)
    Initializes the Amegic object, i.e. in this method this object is instantiated. By calling its
    method InitializeProcesses the processes will be read in and their amplitudes and libraries
    will be constructed and eventually stored in libraries. If the libraries are already
    available, i.e. constructed, compiled and linked, obviuosly this step of storing them is
    unneccessary. In this case, Amegic and te Matrix_Element_Handler are ready for the run
    and 1 is returned. If some libraries still have to be compiled, 0 is returned, resulting in 
    a controlled end of the run.
  */
  /*!
    \fn bool Matrix_Element_Handler::UnweightedEvent()
    Returns .true. if the generation of an unweighted event succeeded in a user-defined number of
    trials.
  */
  /*!
    \fn double Matrix_Element_Handler::WeightedEvent()
    Returns the weight of a new weighted event in \f$1/GeV^2\f$.
  */
  /*!
    \fn Matrix_Element_Handler::Matrix_Element_Handler(std::string,std::string,
                                                       MODEL::Model_Base *,
	    		                               BEAM::Beam_Spectra_Handler *,
                                                       PDF::ISR_Handler *)
    In the constructor the Data_Read object is initialized and the way to calculate
    matrix elements is read in. At the moment this is either through Amegic or through some 
    simple matrix lements of EXTRAXS. These matrix elements are immediately instantiated,
    via either InitializeAmegic or through InitializeSimpleXS. Also, the mode for single 
    event generation is defined to be either "Unweighted" or "Weighted".
  */
  /*!
    \fn Matrix_Element_Handler::~Matrix_Element_Handler()
    Here the initialized matrix elements are deleted.
  */
  /*!
    \fn bool Matrix_Element_Handler::CalculateTotalXSecs()
    This method evaluates the total cross sections for the processes initialized before.
    For the case of Amegic an option is available - but not yet fully implemented and tested -
    to read in the total cross sections, maxima and the a-priori weights for the multi-channel 
    integrator as well as the random number seed. It can not be overstressed that this is quite
    a dangerous option and that it is hte users duty to ensure that the parameters of the
    actual and reference run coincide.
  */
  /*!
    \fn bool Matrix_Element_Handler::PrepareXSecTables()
    Empty at the moment. We have to think whether we need this.
  */
  /*!
    \fn bool Matrix_Element_Handler::LookUpXSec(double,bool,std::string)
    Empty at the moment. We have to think whether we need this.
  */
  /*!
    \fn bool Matrix_Element_Handler::GenerateOneEvent()
    Generates one event, according to m_eventmode this is either a weighted or an unweighted event.
    This is done by the corresponding methods.
  */
  /*!
    \fn bool Matrix_Element_Handler::GenerateSameEvent()
    Generates a new kinematical configuration for a given flavour composition of an event. This is
    important if the matrix elements are merged with the parton showers.
  */
  /*!
    \fn std::string Matrix_Element_Handler::SignalGenerator()
    Returns the user-defined name of the signal process generator he chose.
  */
  /*!
    \fn std::string Matrix_Element_Handler::Name()
    The name of the initialized matrix elements.
  */
  /*!
    \fn std::string Matrix_Element_Handler::ProcessName()
    The name of the actual process for which an event was generated through GenerateOneEvent() or
    GenerateSameEvent(). 
  */
  /*!
    \fn int Matrix_Element_Handler::Nin()
    The number of incoming particles of the actual process for which an event was generated through 
    GenerateOneEvent() or GenerateSameEvent(). 
  */
  /*!
    \fn int Matrix_Element_Handler::Nout()
    The number of outgoing particles of the actual process for which an event was generated through 
    GenerateOneEvent() or GenerateSameEvent(). 
  */
  /*!
    \fn ATOOLS::Vec4D * Matrix_Element_Handler::Momenta()
    The momenta of all particles of the actual process for which an event was generated through 
    GenerateOneEvent() or GenerateSameEvent(). 
  */
  /*!
    \fn ATOOLS::Flavour * Matrix_Element_Handler::Flavs()
    The flavours of all particles of the actual process for which an event was generated through 
    GenerateOneEvent() or GenerateSameEvent(). 
  */
  /*!
    \fn int Matrix_Element_Handler::NumberOfDiagrams()
    This is a method to extract the number of diagrams if the matrix elements generator of choice 
    is Amegic.
  */
  /*!
    \fn AMEGIC::Point * Matrix_Element_Handler::GetDiagram(int)
    This is a method to extract an actual diagrams if the matrix elements generator of choice 
    is Amegic.
  */
  /*!
    \fn EXTRAXS::XS_Base * Matrix_Element_Handler::GetXS()
    This is a method to extract the EXTRAXS::XS_Base in case the simple matrix elements are the user's
    choice.
  */
  /*!
    \fn int Matrix_Element_Handler::MaxJets()
    The maximal number of jets of all hard processes.
  */
  /*!
    \fn double Matrix_Element_Handler::Weight()
    The weight of the last weighted event.
  */    

  /*!
    \fn   PDF::ISR_Handler   * Matrix_Element_Handler::GetISR_Handler()
    This method allows access to the ISR_Handler and the information stored there, e.g.
    which leg encountered ISR, what is the frame of the beams etc.
  */
  /*!
    \fn   int Matrix_Element_Handler::UseSudakovWeight()
    Returns whether to use sudakov weights in single event generation (usually read from ME.dat)
  */

}
#endif

