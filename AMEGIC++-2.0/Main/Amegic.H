#ifndef Amegic_h
#define Amegic_h

#include "All_Processes.H"

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

    int                          m_nmax,m_count;

    All_Processes              * p_procs;
    Interaction_Model_Base     * p_model;
    Topology                   * p_top;

    std::ostream               * p_fifo;

    AORGTOOLS::Data_Read       * p_dataread;
    APHYTOOLS::Selector_Data   * p_seldata;
    BEAM::Beam_Spectra_Handler * p_beam;
    ISR::ISR_Handler           * p_isr;

    void InitializeInteractionModel(MODEL::Model_Base *);
    void ReadInProcessfile(std::string); 
    int  ExtractFlavours(APHYTOOLS::Flavour *&,Pol_Info *&,std::string);
    void Shorten(std::string&);
  public :
    Amegic(std::string,std::string,MODEL::Model_Base * _model);
    ~Amegic();
    bool InitializeProcesses(BEAM::Beam_Spectra_Handler *,ISR::ISR_Handler *);
    bool CalculateTotalXSec();
    bool PrepareXSecTables();
    bool LookUpXSec(double,bool,std::string);
    void SingleEvents();
    void           FifoOutput(double = 1.0) ;
    //!The decay tables of all unstable particles are calculated
    void DecCalc(); 
 };
}

#endif
