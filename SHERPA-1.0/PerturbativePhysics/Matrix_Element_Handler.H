#ifndef Matrix_Element_Handler_h
#define Matrix_Element_Handler_h

#include <set> 
#include <iterator> 
#include <iostream> 

#include "Amegic.H"
#include "Process_Base.H"
#include "SimpleXSecs.H"
#include "XS_Base.H"

#include "Model_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Parton_List.H"
#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Matrix_Element_Handler {
  private :
    std::string                   m_dir,m_file,m_signalgenerator;
    std::string                   m_name;
    std::string                   m_readin;
    bool                          m_eventmode;
    int                           m_mode;

    AORGTOOLS::Data_Read        * p_dataread;
    AMEGIC::Amegic              * p_amegic;
    EXTRAXS::SimpleXSecs        * p_simplexs;

    bool RescaleJetrates() ;

    int InitializeSimpleXS(MODEL::Model_Base *,
			   BEAM::Beam_Spectra_Handler *,ISR::ISR_Handler *);
    int InitializeAmegic(MODEL::Model_Base *,
			 BEAM::Beam_Spectra_Handler *,ISR::ISR_Handler *);
    bool UnweightedEvent();
    bool WeightedEvent();
  public :
    Matrix_Element_Handler(std::string,std::string,MODEL::Model_Base *,
			   BEAM::Beam_Spectra_Handler *,ISR::ISR_Handler *);
    ~Matrix_Element_Handler();

    bool InitializeProcesses(BEAM::Beam_Spectra_Handler *,ISR::ISR_Handler *);
    bool CalculateTotalXSecs();
    bool PrepareXSecTables();
    bool LookUpXSec(double,bool,std::string);

    bool GenerateOneEvent();
    bool GenerateSameEvent();

    std::string          SignalGenerator();
    std::string          ProcessName();
    std::string          Name();
    int                  Nin();
    int                  Nout();
    int                  MaxJets();
    AMATOOLS::Vec4D    * Momenta();
    APHYTOOLS::Flavour * Flavs();

    int                  NumberOfDiagrams();
    AMEGIC::Point      * GetDiagram(int);
    EXTRAXS::XS_Base   * GetXS();
  };
}
#endif

