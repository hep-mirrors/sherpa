#ifndef Amisic_H
#define Amisic_H

#include "Grid_Handler.H"

#include "SimpleXSecs.H"
#include "QCD_Processes.H"
#include "Blob_List.H"

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
    typedef Grid_Creator_Base<GridArgumentType,GridResultType> GridCreatorBaseType;
    typedef Grid_Creator<GridArgumentType,GridResultType> GridCreatorType;

  private:

    std::vector<ATOOLS::Blob_List> m_blobs;
    std::vector<std::string>       m_filename, m_processname;
    std::vector<bool>              m_create;
    std::vector<GridHandlerType*>  m_differential;

    GridHandlerType *m_total;
    
    std::string m_inputdirectory, m_inputfile, m_environmentfile, m_xsfile;
    std::string m_outputdirectory;

    ATOOLS::Blob *CreateBlob(ATOOLS::Flavour flavour[4]);

    void ClearVector();
    bool HaveBlob(ATOOLS::Blob *blob);
    void FillMode(EXTRAXS::QCD_Processes::Mode mode);

    bool ReadInData();
    bool CreateGrid(ATOOLS::Blob_List& bloblist,
		    std::string& filename,std::string& processname);
    bool CalculateTotal();

  public:

    // constructor
    Amisic();

    // destructor
    ~Amisic();

    // setters 
    void SetInputFile(std::string _m_inputfile);
    void SetEnvironmentFile(std::string _m_environmentfile);
    void SetInputDirectory(std::string _m_inputdirectory);
    void SetOutputDirectory(std::string _m_outputdirectory);

    // getters
    std::string GetInputFile();
    std::string GetEnvironmentFile();
    std::string GetInputDirectory();
    std::string GetOutputDirectory();

    // member functions
    bool Initialize(std::string tempidir=ATOOLS::nullstring,
		    std::string tempifile=ATOOLS::nullstring,bool creategrid=true);
    // ATOOLS::Blob_List CreateProcesses();

  }; // end of class Amisic

} // end of namespace AMISIC

#include "Amisic.IC"

#endif
