#ifndef Amisic_H
#define Amisic_H

#include "Grid_Handler.H"

#include "SimpleXSecs.H"
#include "Environment.H"
#include "QCD_Processes_C.H"
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

      typedef std::vector<GridArgumentType> GridArgumentVector;
      typedef std::vector<GridResultType> GridResultVector;

      typedef Grid_Handler<GridArgumentType,GridResultType> GridHandlerType;
      typedef std::vector<GridHandlerType*> GridHandlerVector;
    
      typedef ATOOLS::Data_To_Function<GridArgumentType,GridResultType> GridFunctionType;

    private:
      
      GridHandlerVector     p_gridhandler;
      EXTRAXS::SimpleXSecs *p_processes;

      int m_criterion;
      bool m_storemax;

      ATOOLS::Mom_Data m_initialdata;
      GridResultType m_lastxs, m_lastmax;

      ATOOLS::Axis<GridArgumentType> *p_xaxis;
      ATOOLS::Axis<GridResultType> *p_yaxis;

      std::string m_xsextension, m_maxextension;

    public:

      // constructor
      Grid_Creator(GridHandlerVector _p_gridhandler,EXTRAXS::SimpleXSecs *_p_processes);
      
      // setters
      void SetXSExtension(std::string _m_xsextension);
      void SetMaxExtension(std::string _m_maxextension);

      // getters
      std::string XSExtension();
      std::string MaxExtension();

      GridHandlerType *RelevantGrid();
      
      // member functions
      bool InitializeCalculation();
      unsigned int CreateSinglePoint(GridArgumentType *boundary,
				     bool newpoint,bool force);

      bool CreateInitialGrid();
      bool CreateOptimizedGrid();

      bool ReadInArguments(std::string tempifile=ATOOLS::nullstring,
			   std::string tempipath=ATOOLS::nullstring);
      bool WriteOutGrid(std::vector<std::string> addcomments=std::vector<std::string>(0),
			std::string tempofile=ATOOLS::nullstring,
			std::string tempopath=ATOOLS::nullstring);

    }; // end of class Amisic::Grid_Creator

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

    long int m_particlecounter;

    GridArgumentType m_start, m_stop;
    
    std::string m_inputpath, m_inputfile, m_environmentfile, m_xsfile;
    std::string m_outputpath;
    std::string m_xsextension, m_maxextension;

    AMEGIC::Environment  *p_environment;
    EXTRAXS::SimpleXSecs *p_processes;
    FSRChannel           *p_fsrinterface;

    ATOOLS::Blob *GetBlob(ATOOLS::Flavour flavour[4]);

    void CleanUp();
    void FillMode(EXTRAXS::QCD_Processes_C::Mode mode);
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

    void SetInputPath(std::string _m_inputdirectory);
    void SetInputFile(std::string _m_inputfile);
    void SetEnvironmentFile(std::string _m_environmentfile);
    void SetOutputPath(std::string _m_outputdirectory);

    // getters
    GridArgumentType Start();
    GridArgumentType Stop();

    std::string InputPath();
    std::string InputFile();
    std::string EnvironmentFile();
    std::string OutputPath();

    // member functions
    bool Initialize(std::string tempidir=ATOOLS::nullstring,
		    std::string tempifile=ATOOLS::nullstring,bool creategrid=true);
    ATOOLS::Blob_List *CreateProcesses();

  }; // end of class Amisic

} // end of namespace AMISIC

#include "Amisic.IC"

#endif
