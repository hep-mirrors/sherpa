#ifndef Amegic_h
#define Amegic_h

#include "All_Processes.H"
#include "Selector.H"
#include "Flavour.H"
#include "MyString.H"

#include "Beam_Handler.H"
#include "ISR_Handler.H"


namespace AMEGIC {

  class Amegic {
    int    gen_str;    // 0 = no strings, 1 = strings, 2 = strings + libs
    std::string                path;
    std::string                resdir;

    bool                       beamon, isron;
    int                      * beamtypes, * isrtypes;
    double                   * splimits;
    int                        arrows;
    APHYTOOLS::Flavour       * bunches, * beams, * partons;
    Pol_Info                 * plpartons;
    double                   * polbunches, * polbeams;

    APHYTOOLS::Selector_Data * seldata;
    ISR::ISR_Handler         * isr;
    BEAM::Beam_Handler       * beam;

    All_Processes            * procs;

    Topology                 * top;

    int runmode;
    
      //int ExtractFlavours(APHYTOOLS::Flavour*&,std::string);
    int ExtractFlavours(APHYTOOLS::Flavour*&,Pol_Info*&,std::string);
    int ExtractFlavours(APHYTOOLS::Flavour*&,double*&,std::string);
    int ReadProcesses(std::string);
  public:
    Amegic(std::string path,ISR::ISR_Handler * _isr=0,BEAM::Beam_Handler * _beam = 0);
    ~Amegic();
    //! Calculates the total cross section of all involved processes.
    bool           CalculateTotalXSec();
    bool           PrepareXSecTables();
    bool           LookUpXSec(double,bool,std::string);
    void           SingleEvents();
    bool           InitializeProcesses(int _mode=AMPLITUDE_MODE);
    Process_Base * Processes();

    void           SetResDir(std::string _resdir) { resdir = _resdir; }
    std::string    ResDir()                       { return resdir; }
  };
  
  inline Process_Base * Amegic::Processes() { return procs; };
}

#endif

//! The main steering class.(ready)
/*!
  It handles and calculates all processes to be evaluated or
  taken into account.
*/
    //! The list of Single_Process(es)/Process_Group(s).
    //! All possible topologies.
    //! Extracts the makronames for the processes or the flavours of specific ones.
    /*!
      The raw text given in the Processes.dat is translated into the numbers
      and Flavours of the incoming and outgoing particles.

      I have defined the following makros :  
      IS                              FS               Description
      {11,-11},{13,-13},{21,21}       2_3_...n_jets    strong jets up to n
    */ 
    //! Link to predefined process classes
    /*
      int ProcessClasses(std::string,APHYTOOLS::Flavour *&,ISR::ISR_Handler *,BEAM::Beam_Handler *);
      
      int InitMultiJets(std::string,int,APHYTOOLS::Flavour *&,
      ISR::ISR_Handler *,BEAM::Beam_Handler *,std::string);
      int InitFourFermions(std::string,int,APHYTOOLS::Flavour *&,
      ISR::ISR_Handler *,BEAM::Beam_Handler *,std::string);
      int InitDY(std::string,int,APHYTOOLS::Flavour *&,
      ISR::ISR_Handler *,BEAM::Beam_Handler *,std::string);
    */
    //! Reads Processes from Processes.dat.
    /*!
      Reads in the list of processes from Process.dat, extracts the Flavours 
      with the help of ExtractFlavours, and groups the processes into 
      All_Processes.For specific keywords full process-classes are designed.
      They are filled via ProcessClasses.
      \sa ExtractFlavours(), ProcessClasses()
    */ 
    //! The constructor of this class
    /*!
      Runs the whole initialization procedure of the program following these steps:
      - Set the Model file
      - Initialize the Model
      - Set all vertices within the Model
      - Read in the processes
      - Initialize the different topologies
      - Calculate the decay width and the branching fractions of the particles using the
        Decay_Handler. 
      - Initializes the Process	 
      \sa ReadProcesses(), Vertex, Topology
    */ 
