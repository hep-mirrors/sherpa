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

    ATOOLS::Data_Read       * p_dataread;
    ATOOLS::Selector_Data   * p_seldata;
    BEAM::Beam_Spectra_Handler * p_beam;
    PDF::ISR_Handler           * p_isr;

    void InitializeInteractionModel(MODEL::Model_Base *);
    void ReadInProcessfile(std::string); 
    int  ExtractFlavours(ATOOLS::Flavour *&,Pol_Info *&,std::string);
    void Shorten(std::string&);
  public :
    Amegic(std::string,std::string,MODEL::Model_Base * _model);
    ~Amegic();
    bool        InitializeProcesses(BEAM::Beam_Spectra_Handler *,PDF::ISR_Handler *);
    bool        InitializeDecays();
    bool        CalculateTotalXSec(std::string _resdir=std::string(""));
    bool        CalculateBranchingWidths(std::string _resdir=std::string(""));
    void        SetResDir(std::string);

    bool        UnweightedEvent();
    double      WeightedEvent();
    bool        SameEvent();
    double      SameWeightedEvent();
    Process_Base       * GetProcess()       { return p_procs->Selected(); }
    Process_Base * Processes();
    std::string          ProcessName()      { return p_procs->Selected()->Name(); }
    int                  Nin()              { return p_procs->Selected()->Nin(); }
    int                  Nout()             { return p_procs->Selected()->Nout(); }
    ATOOLS::Flavour * Flavs()            { return p_procs->Selected()->Flavs(); }
    ATOOLS::Vec4D    * Momenta()          { return p_procs->Selected()->Momenta(); }
    int                  NumberOfDiagrams() { return p_procs->Selected()->NumberOfDiagrams(); }
    AMEGIC::Point      * Diagram(int _diag) { return p_procs->Selected()->Diagram(_diag); }
    int                  MaxJets()          { return m_maxjet; }


    void SingleEvents();
    void FifoOutput(double = 1.0) ;

    //!The decay tables of all unstable particles are calculated
    void DecCalc(); 

    // not used right now !!!
    bool        LookUpXSec(double,bool,std::string);
    bool        PrepareXSecTables();
 };

  inline Process_Base * Amegic::Processes() { return p_procs; };

}

#endif


