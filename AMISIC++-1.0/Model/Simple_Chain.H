#ifndef Simple_Chain_H
#define Simple_Chain_H

#include "MI_Base.H"
#include "Grid_Creator.H"
#include "Data_To_Function.H"
#include "Profile_Function_Base.H"
#include "Fixed_Variable_Channel.H"
#include "Simple_XS.H"
#include "Environment.H"
#ifdef USING__Sherpa
#include "Remnant_Base.H"
#endif

namespace SHERPA { class Remnant_Base; }

namespace AMISIC {

  class Simple_Chain: public MI_Base {
  public:

    typedef ATOOLS::Data_To_Function<double,double> Grid_Function_Type;

    typedef Amisic_Histogram<double> Amisic_Histogram_Type;
    typedef std::map<std::string,Amisic_Histogram_Type*>  Amisic_Histogram_Map;

    typedef PHASIC::Fixed_Variable_Channel FSR_Channel;

    typedef std::map<std::string,EXTRAXS::XS_Base*> Process_Map;
    typedef std::multimap<double,std::string,
			  std::greater<double> >    Sort_Map;

  private:

    Amisic_Histogram_Map   m_differentials;
    Grid_Function_Type    *p_total;
    Amisic_Histogram_Type *p_differential;
    Grid_Creator          *p_gridcreator;

    Process_Map m_processmap;

    double m_norm, m_enhance, m_maxreduction;

    std::string m_environmentfile, m_xsfile;
    std::string m_xsextension, m_mcextension;

    EXTRAXS::Simple_XS *p_processes;
    FSR_Channel        *p_fsrinterface;

    AMISIC::Environment        *p_environment;
    MODEL::Model_Base          *p_model;
    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;

    Profile_Function_Base *p_profile;

    std::string m_selected;

    unsigned int m_nflavour, m_maxtrials, m_check;
    unsigned int m_scalescheme, m_kfactorscheme, m_vegas, m_foam;
    double       m_sigmahard, m_error, m_ecms;
    bool         m_external, m_regulate;

    std::string         m_regulator; 
    std::vector<double> m_regulation;

    SHERPA::Remnant_Base *p_remnants[2];

    void CleanUp();

    EXTRAXS::XS_Group *FindPDFGroup(const size_t nin,const size_t nout,
				    const ATOOLS::Flavour *flavours);

    bool GeneratePathName();
    void OrderFlavours(ATOOLS::Flavour *flavours);
    bool AddProcess(EXTRAXS::XS_Group *const group,
		    const ATOOLS::Flavour *flavours);
    bool ReadInData();

    bool CreateGrid();
    bool InitializeProcessList();
    void CalculateSigmaND();
    bool CalculateTotal();
    bool SetUpInterface();

    bool CheckConsistency(EXTRAXS::XS_Group *const group,
			  Grid_Function_Type *const grid,
			  const double min,const double max,
			  const double integral);

    bool FillBlob(ATOOLS::Blob *blob);

    void PrepareTerminate();

  public:

    // constructor
    Simple_Chain();
    Simple_Chain(MODEL::Model_Base *_p_model,
		 BEAM::Beam_Spectra_Handler *_p_beam,
		 PDF::ISR_Handler *_p_isr);

    // destructor
    ~Simple_Chain();

    // member functions
    bool Initialize();
    void Update(const MI_Base *mibase);
    void Reset();

    bool DiceEnhanceFactor();
    bool DiceOrderingParameter();

    bool VetoProcess(ATOOLS::Blob *blob);
    bool DiceProcess();

    // inline functions
    inline double Norm() { return m_norm; }

    inline std::string EnvironmentFile() { return m_environmentfile; }

    inline Profile_Function_Base *ProfileFunction() { return p_profile; }

    inline void SetNorm(double norm) { m_norm=norm; }

    inline void SetEnvironmentFile(std::string environmentfile)
    { m_environmentfile=environmentfile; }

  }; // end of class Simple_Chain

} // end of namespace AMISIC

#endif
