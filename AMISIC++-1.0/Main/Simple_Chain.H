#ifndef Simple_Chain_H
#define Simple_Chain_H

#include "MI_Base.H"
#include "Grid_Creator.H"
#include "QCD_Processes.H"

namespace AMISIC {
  
  class Simple_Chain: public MI_Base {
  public:

    typedef double GridArgumentType;
    typedef double GridResultType;

    typedef Grid_Handler<GridArgumentType,GridResultType> GridHandlerType;
    typedef Grid_Creator<GridArgumentType,GridResultType> GridCreatorType;
    typedef Grid_Creator_Base<GridArgumentType,GridResultType> GridCreatorBaseType;
    typedef Data_To_Function<GridArgumentType,GridResultType> GridFunctionType;
    typedef std::vector<GridHandlerType*> GridHandlerVector;

    typedef PHASIC::Fixed_Variable_Channel<GridArgumentType> FSRChannel;

  private:

    std::vector<ATOOLS::Blob_List> m_blobs;
    std::vector<std::string>       m_filename, m_comments;

    std::vector<GridFunctionType*>          m_differential;
    std::map<std::string,GridFunctionType*> m_maxima;

    GridResultType    m_norm;
    GridFunctionType *p_total;

    GridCreatorType *p_gridcreator;

    std::string m_environmentfile, m_xsfile;
    std::string m_xsextension, m_maxextension, m_mcextension;

    EXTRAXS::Simple_XS *p_processes;
    FSRChannel         *p_fsrinterface;

    AMEGIC::Environment        *p_environment;
    MODEL::Model_Base          *p_model;
    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;

    unsigned int m_selected, m_nflavour, m_scalescheme;
    unsigned int m_kfactorscheme, m_maxtrials;
    double       m_sigmahard, m_error, m_ecms;
    bool         m_external;

    ATOOLS::Blob *GetBlob(ATOOLS::Flavour *flavour);

    void CleanUp();
    bool HaveBlob(ATOOLS::Blob *blob);
    void FillMode(EXTRAXS::QCD_Processes::Mode mode);
    void SetMaximum(EXTRAXS::XS_Base *process); 

    bool ReadInData();
    bool CreateGrid(ATOOLS::Blob_List& bloblist,std::string& filename);
    bool CalculateTotal();
    void CalculateSigmaND();
    bool InitializeBlobList();

    bool FillBlob(ATOOLS::Blob *blob);
    void PrepareTerminate();

  public:

    // constructor
    Simple_Chain();
    Simple_Chain(MODEL::Model_Base *_p_model,
		 BEAM::Beam_Spectra_Handler *_p_beam,PDF::ISR_Handler *_p_isr);

    // destructor
    ~Simple_Chain();

    // setter
    void SetEnvironmentFile(std::string _m_environmentfile);
    void SetNorm(GridResultType _m_norm);

    // getter
    std::string    EnvironmentFile();
    GridResultType Norm();

    // member functions
    bool Initialize();
    void Update(const MI_Base *mibase);
    void Reset();

    bool DiceOrderingParameter();
    bool DiceProcess();

  }; // end of class Simple_Chain

} // end of namespace AMISIC

#include "Simple_Chain.inl.H"

#endif
