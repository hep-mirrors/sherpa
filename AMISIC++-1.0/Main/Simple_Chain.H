#ifndef Simple_Chain_H
#define Simple_Chain_H

#include "MI_Base.H"
#include "Grid_Creator.H"

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

    std::vector<ATOOLS::Blob_List>  m_blobs;
    std::vector<std::string>        m_filename, m_processname;
    std::vector<bool>               m_create;
    std::vector<GridFunctionType*>  m_differential, m_maximum;
    GridFunctionType               *p_total;

    std::string m_environmentfile, m_xsfile;
    std::string m_xsextension, m_maxextension;

    EXTRAXS::SimpleXSecs *p_processes;
    FSRChannel           *p_fsrinterface;

    AMEGIC::Environment        *p_environment;
    MODEL::Model_Base          *p_model;
    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;

    bool m_external;

    unsigned int m_selected;

    ATOOLS::Blob *GetBlob(ATOOLS::Flavour *flavour);

    void CleanUp();
    bool HaveBlob(ATOOLS::Blob *blob);
    void FillMode(EXTRAXS::QCD_Processes_C::Mode mode);

    bool ReadInData();
    bool CreateGrid(ATOOLS::Blob_List& bloblist,
		    std::string& filename,std::string& processname);
    bool CalculateTotal();
    bool InitializeBlobList();

    bool FillBlob(ATOOLS::Blob *blob);

  public:

    // constructor
    Simple_Chain();
    Simple_Chain(MODEL::Model_Base *_p_model,
		 BEAM::Beam_Spectra_Handler *_p_beam,PDF::ISR_Handler *_p_isr);

    // destructor
    ~Simple_Chain();

    // setter
    void SetEnvironmentFile(std::string _m_environmentfile);

    // getter
    std::string EnvironmentFile();

    // member functions
    bool Initialize();
    void Update(MI_Base *mibase);
    void Reset();

    bool DiceOrderingParameter();
    bool DiceProcess();

  }; // end of class Simple_Chain

} // end of namespace AMISIC

#include "Simple_Chain.IC"

#endif
