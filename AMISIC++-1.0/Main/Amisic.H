#ifndef Amisic_H
#define Amisic_H

#include "Grid_Handler.H"

#include "SimpleXSecs.H"
#include "Environment.H"
#include "QCD_Processes.H"
#include "Blob_List.H"
#include "Fixed_Variable_Channel.H"

namespace AMISIC {
  
  class Amisic {
  public:

    template <class Argument_Type,class Result_Type>
    class Grid_Creator: public AMISIC::Grid_Creator_Base<Argument_Type,Result_Type> {
    public:
      
      typedef Argument_Type GridArgumentType;
      typedef Result_Type GridResultType;

      typedef Grid_Handler<GridArgumentType,GridResultType> GridHandlerType;

    private:
      
      EXTRAXS::SimpleXSecs *p_processes;

      int m_criterion;
      ATOOLS::Mom_Data m_initialdata;
      ATOOLS::Axis<GridArgumentType> *p_xaxis;

    public:

      // constructor
      Grid_Creator(GridHandlerType *_p_gridhandler,
		   EXTRAXS::SimpleXSecs *_p_processes);
      
      // member functions
      bool InitializeCalculation();
      GridResultType CalculateSingleValue(GridArgumentType nextleft,
					  GridArgumentType nextright);
      
    }; // end of class Amisic::Grid_Creator

    typedef double GridArgumentType;
    typedef double GridResultType;

    typedef Grid_Handler<GridArgumentType,GridResultType> GridHandlerType;
    typedef Grid_Creator<GridArgumentType,GridResultType> GridCreatorType;
    typedef Grid_Creator_Base<GridArgumentType,GridResultType> GridCreatorBaseType;
    typedef Data_To_Function<GridArgumentType,GridResultType> GridFunctionType;
    
    typedef PHASIC::Fixed_Variable_Channel<GridArgumentType> FSRChannel;

  private:

    std::vector<ATOOLS::Blob_List>  m_blobs;
    std::vector<std::string>        m_filename, m_processname;
    std::vector<bool>               m_create;
    std::vector<GridFunctionType*>  m_differential;
    GridFunctionType               *p_total;

    GridArgumentType m_start, m_stop;
    
    std::string m_inputdirectory, m_inputfile, m_environmentfile, m_xsfile;
    std::string m_outputdirectory;

    AMEGIC::Environment  *p_environment;
    EXTRAXS::SimpleXSecs *p_processes;
    FSRChannel           *p_fsrinterface;

    ATOOLS::Blob *GetBlob(ATOOLS::Flavour flavour[4]);

    void CleanUp();
    void FillMode(EXTRAXS::QCD_Processes::Mode mode);
    bool HaveBlob(ATOOLS::Blob *blob);

    bool ReadInData();
    bool CreateGrid(ATOOLS::Blob_List& bloblist,
		    std::string& filename,std::string& processname);
    bool CalculateTotal();
    bool InitializeBlobList();

    ATOOLS::Blob     *CreateBlob(unsigned int bloblist,double ran);
    ATOOLS::Blob     *DiceProcess(GridArgumentType parameter,double ran[2]);
    GridArgumentType  DiceParameter(GridArgumentType &last,double ran);

  public:

    // constructor
    Amisic();

    // destructor
    ~Amisic();

    // setters 
    void SetStart(GridArgumentType _m_start);
    void SetStop(GridArgumentType _m_stop);

    void SetInputFile(std::string _m_inputfile);
    void SetEnvironmentFile(std::string _m_environmentfile);
    void SetInputDirectory(std::string _m_inputdirectory);
    void SetOutputDirectory(std::string _m_outputdirectory);

    // getters
    GridArgumentType Start();
    GridArgumentType Stop();

    std::string InputFile();
    std::string EnvironmentFile();
    std::string InputDirectory();
    std::string OutputDirectory();

    // member functions
    bool Initialize(std::string tempidir=ATOOLS::nullstring,
		    std::string tempifile=ATOOLS::nullstring,bool creategrid=true);
    ATOOLS::Blob_List *CreateProcesses();

  }; // end of class Amisic

} // end of namespace AMISIC

#include "Amisic.IC"

#endif
