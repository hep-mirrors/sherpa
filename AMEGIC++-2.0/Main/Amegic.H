#ifndef Amegic_h
#define Amegic_h

#include "All_Processes.H"
#include "All_Decays.H"

#include "Model_Base.H"
#include "Interaction_Model_Base.H"

#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Selector.H"

#include "Flavour.H"
#include "Pol_Info.H"
#include "Data_Read.H"


namespace AMEGIC {
  class Amegic {
  private :
    std::string                  m_path,m_file;
    std::string                  m_respath,m_resfile;

    int                          m_nmax,m_count,m_maxjet;

    All_Processes              * p_procs;
    All_Decays                 * p_decs;
    Interaction_Model_Base     * p_model;
    Topology                   * p_top;

    std::ostream               * p_fifo;

    ATOOLS::Data_Read          * p_dataread;
    ATOOLS::Selector_Data      * p_seldata;
    BEAM::Beam_Spectra_Handler * p_beam;
    PDF::ISR_Handler           * p_isr;

    void InitializeInteractionModel(MODEL::Model_Base *);
    void ReadInProcessfile(std::string); 
    int  ExtractFlavours(ATOOLS::Flavour *&,Pol_Info *&,std::string);
    void Shorten(std::string&);
  public :
    Amegic(std::string,std::string,MODEL::Model_Base *);
    ~Amegic();
    bool        InitializeProcesses(BEAM::Beam_Spectra_Handler *,PDF::ISR_Handler *);
    bool        InitializeDecays(bool);
    bool        CalculateTotalXSec(std::string _resdir=std::string(""), int mode=0);
    bool        CalculateBranchingWidths(std::string _resdir=std::string(""));
    void        SetResDir(std::string);

    bool              UnweightedEvent();
    ATOOLS::Blob_Data_Base * WeightedEvent();
    bool              SameEvent();
    ATOOLS::Blob_Data_Base * SameWeightedEvent();
    PHASIC::Integrable_Base  *const GetProcess() const { return p_procs->Selected(); }
    Process_Base    * Processes();
    All_Decays      * GetAllDecays()     { return p_decs; }
    std::string       ProcessName()      { return p_procs->Selected()->Name(); }
    int               NIn()              { return p_procs->Selected()->NIn(); }
    int               NOut()             { return p_procs->Selected()->NOut(); }
    const ATOOLS::Flavour * Flavours()   { return p_procs->Selected()->Flavours(); }
    const ATOOLS::Vec4D   * Momenta()    { return p_procs->Selected()->Momenta(); }
    int               InSwaped()         { return p_procs->Selected()->InSwaped(); }
    int               OrderStrong()      { return static_cast<Process_Base*>(p_procs->Selected())->OrderStrong(); }
    int               OrderEWeak()       { return static_cast<Process_Base*>(p_procs->Selected())->OrderEWeak(); }
    int               NumberOfDiagrams() { return static_cast<Process_Base*>(p_procs->Selected())->NumberOfDiagrams(); }
    AMEGIC::Point   * Diagram(int _diag) { return static_cast<Process_Base*>(p_procs->Selected())->Diagram(_diag); }
    int               MaxJets()          { return m_maxjet; }


    void SingleEvents();
    void FifoOutput(double = 1.0) ;
    double OverflowStatistics(Process_Base * proc=NULL,int level=0);

    //!The decay tables of all unstable particles are calculated
    void DecCalc(); 

    // not used right now !!!
    bool        LookUpXSec(double,bool,std::string);
    bool        PrepareXSecTables();
    
    void AddEvent(const double xs,const double validxs,const int ncounts);
    
 };

  inline Process_Base * Amegic::Processes() { return p_procs; }

}

#endif


