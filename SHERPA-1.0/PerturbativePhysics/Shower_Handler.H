#ifndef Shower_Handler_h
#define Shower_Handler_h 

#include "Model_Base.H"
#include "ISR_Handler.H"
#include "Adicic.H"
#include "Hard_Interface.H"
#include "Data_Read.H"

#include <string>
#include <iostream> 

namespace APACIC { class Tree; }
namespace ADICIC { class Chain; }

namespace SHERPA {
  class Shower_Handler {
  private :
    std::string                  m_dir,m_file;
    int                          m_maxjetnumber;


    std::string                  m_showergenerator;
    bool                         m_isrshowerswitch,m_fsrshowerswitch;
    ATOOLS::Data_Read          * p_dataread;
    APACIC::Hard_Interface     * p_apacic;
    ADICIC::Adicic             * p_adicic;
    PDF::ISR_Handler           * p_isr_handler;       
  public :
    Shower_Handler(std::string,std::string,MODEL::Model_Base *,
		   PDF::ISR_Handler *,int);
    ~Shower_Handler();
    int             PerformShowers(int,double,double);
    int             PerformDecayShowers(bool);
    void            FillBlobs(ATOOLS::Blob_List *);
    void            FillDecayBlobs(ATOOLS::Blob_List *);
    void            CleanUp();
    std::string     ShowerGenerator() { return m_showergenerator; }
    bool            ISROn()           { return m_isrshowerswitch; }
    bool            FSROn()           { return m_fsrshowerswitch; }
    unsigned int    MaxJetNumber()    { return m_maxjetnumber;    }
    PDF::ISR_Handler * GetISRHandler() { return p_isr_handler; }
    APACIC::Tree  * GetFinTree(); 
    APACIC::Tree ** GetIniTrees();
    APACIC::Hard_Interface * GetApacic() { return p_apacic; }
    ADICIC::Chain * GetChain();
    void SetJetvetoPt2(const double pt2) { p_apacic->SetJetvetoPt2(pt2); }
    void SetFactorisationScale(const double scale) { p_apacic->SetFactorisationScale(scale); }
  };


  /*!
    \class Shower_Handler
    \brief A wrapper around shower generators. 

    This class provides a generic interface of shower generators for multiple
    emission of secondary partons during the jet evolution phase with the rest
    of the SHERPA framework. 

    At the moment, only one shower generator is available, namely Apacic. Its shower 
    treatment is quite similar to the one employed by Pythia in the sense that the 
    same order parameter is used for the simulation of multiple parton emission, namely 
    the virtual mass of the decaying parton. To model the effect of colour coherence
    angular ordering is imposed by hand in this framework.

    It can be anticipated that another shower option will be made available in the future.
    However, if users want to interface their own parton shower treatment, they would have
    to modify the shower handler accordingly and provide their version of the interface
    to the matrix elements.
  */
  /*!
    \var std::string Shower_Handler::m_dir
    The directory for the file with the shower data.
  */
  /*!
    \var std::string Shower_Handler::m_file
    The name of the data file with the shower specificiations.
  */
  /*!
    \var std::string Shower_Handler::m_showergenerator
    The name of the shower generator. Right now, the only option is APACIC.
  */
  /*!
    \var int Shower_Handler::m_maxjetnumber
    The maximal number of jets stemming from the Signal_Processes. This is needed
    for the set-up of the interface between matrix elements and the shower.
  */
  /*!
    \var bool Shower_Handler::m_isrshowerswitch
    .true. if the initial state shower is employed.
  */
  /*!
    \var bool Shower_Handler::m_fsrshowerswitch
    .true. if the final state shower is employed.
  */
  /*!
    \var ATOOLS::Data_Read * Shower_Handler::p_dataread
    A pointer to the Data_Read object. Will be passed into the shower generators
    for further reading of the input file.
  */
  /*!
    \var APACIC::Hard_Interface * Shower_Handler::p_apacic
    A pointer to the Apacic object. THis is the only shower handler that is actually
    implemented.
  */
  /*!
    \fn Shower_Handler::Shower_Handler(std::string,std::string,MODEL::Model_Base *,
		                       PDF::ISR_Handler *,int)
    The constructor reads in some basic information concerning the shower and then
    initializes the shower generator of choice. The Data_Read object is passed as an
    argument to the constructor of the shower generator for more detailled specification
    of the shpwer parameters.
  */
  /*!
    \fn Shower_Handler::~Shower_Handler()
    Deletes the shower generator.
  */
  /*!
    \fn int Shower_Handler::PerformShowers(bool)
    Provides a link to the corresponding methods of the shower generators.
  */
  /*!
    \fn void Shower_Handler::FillBlobs(ATOOLS::Blob_List *)
    Blobs for each jet are produced and filled by suitable methods of the shower generator.
  */
  /*!
    \fn void Shower_Handler::CleanUp()
    The binary trees of Apacic are prepared for a new event, i.e. they are cleaned to be filled
    again during the next jet evolution.
  */
  /*!
    \fn APACIC::Tree * Shower_Handler::GetFinTree() 
    This is meant specifically, if Apacic is the shower generator of choice and returns the binary tree
    through which the final state parton shower is organized.
  */
  /*!
    \fn APACIC::Tree ** Shower_Handler::GetIniTrees()
    This is meant specifically, if Apacic is the shower generator of choice and returns the two binary 
    trees through which the initial state parton shower is organized.
  */
  /*!
    \fn std::string Shower_Handler::ShowerGenerator()
    Returns the name of the actual shower generator.
  */
  /*!
    \fn bool Shower_Handler::ISROn()
    Returns .true. is the initial state shower is used.
  */
  /*!
    \fn bool Shower_Handler::FSROn()
    Returns .true. is the final state shower is used.
  */
  /*!
    \fn int Shower_Handler::MaxJetNumber()
    Returns the maximal number of jets stemming from the matrix elements.
  */
}


#endif
