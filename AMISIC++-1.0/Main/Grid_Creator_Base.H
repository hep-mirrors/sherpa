#ifndef Grid_Creator_Base_H
#define Grid_Creator_Base_H

#include "File_IO_Base.H"
#include "Grid_Handler.H"

namespace AMISIC { template <class Argument_Type,class Result_Type>
		   class Grid_Handler;                              }

namespace AMISIC {

  template <class Argument_Type,class Result_Type>
  class Grid_Creator_Arguments {
  public:
    
    typedef Argument_Type GridArgumentType;
    typedef Result_Type GridResultType;
    
  protected:
    
    GridArgumentType m_gridxmin, m_gridxmax, m_griddeltaxmax, m_griddeltaxmin;
    GridResultType   m_gridymin, m_gridymax, m_griddeltaymax, m_griddeltaymin;
    
  public:
    
    std::string m_gridxvariable, m_gridxscaling;
    std::string m_gridyvariable, m_gridyscaling;
    
    bool m_useymax, m_useymin, m_optimize;
    
    bool m_relativexmin, m_relativedeltaxmin;
    bool m_relativeymin, m_relativedeltaymin;
    
    int m_outputlevel;
    
  public:
    
    // constructor
    Grid_Creator_Arguments();
    
    // setters
    void SetGridXMin(const GridArgumentType _m_gridxmin);
    void SetGridXMax(const GridArgumentType _m_gridxmax);
    void SetGridDeltaXMin(const GridArgumentType _m_griddeltaxmin);
    void SetGridDeltaXMax(const GridArgumentType _m_griddeltaxmax);
    void SetGridYMin(const GridResultType _m_gridymin);
    void SetGridYMax(const GridResultType _m_gridymax);
    void SetGridDeltaYMin(const GridResultType _m_griddeltaymin);
    void SetGridDeltaYMax(const GridResultType _m_griddeltaymax);
    
    // getters
    const GridArgumentType GridXMin() const;
    const GridArgumentType GridXMax() const;
    const GridArgumentType GridDeltaXMin() const;
    const GridArgumentType GridDeltaXMax() const;
    const GridResultType   GridYMin() const;
    const GridResultType   GridYMax() const;
    const GridResultType   GridDeltaYMin() const;
    const GridResultType   GridDeltaYMax() const;
    
  }; // end of class Arguments

  template <class Argument_Type,class Result_Type>
  class Grid_Creator_Base: public File_IO_Base,
			   public Grid_Creator_Arguments<Argument_Type,Result_Type> {
  public:
    
    typedef Argument_Type GridArgumentType;
    typedef Result_Type GridResultType;

    typedef std::vector<GridArgumentType> GridArgumentVector;
    typedef std::vector<GridResultType> GridResultVector;

    typedef Grid_Handler<GridArgumentType,GridResultType> GridHandlerType;
    typedef Grid_Creator_Arguments<GridArgumentType,GridResultType> ArgumentsType;
    typedef ATOOLS::Data_To_Function<GridArgumentType,GridResultType> GridFunctionType;

    typedef std::vector<GridHandlerType*> GridHandlerVector;

  protected:

    virtual bool InitializeCalculation();
    virtual unsigned int CreateSinglePoint(GridArgumentType *boundary,
					   bool newpoint=true,bool force=false);

    bool InitializeSingleGrid(GridHandlerType *grid);
    bool OptimizeSingleGrid(GridHandlerType *grid);

    virtual bool CreateInitialGrid();
    virtual bool CreateOptimizedGrid();

    bool ReadSingleArguments(GridHandlerType *grid);
#ifndef DEBUG__Simple_Chain
    bool WriteSingleGrid(GridHandlerType *grid,std::vector<std::string> addcomments);
#endif

  public:

    // constructors
    Grid_Creator_Base();

    // destructor
    virtual ~Grid_Creator_Base();
    
    // setters
    bool SetArguments(GridHandlerType *grid,ArgumentsType _m_arguments);

    // getters
    ArgumentsType Arguments();

    // member functions
    virtual bool ReadInArguments(std::string tempifile=ATOOLS::nullstring,
				 std::string tempipath=ATOOLS::nullstring);
    virtual bool WriteOutGrid(std::vector<std::string> addcomments=std::vector<std::string>(0),
			      std::string tempofile=ATOOLS::nullstring,
			      std::string tempopath=ATOOLS::nullstring);
    
    bool CreateGrid();

#ifdef DEBUG__Simple_Chain
    bool WriteSingleGrid(GridHandlerType *grid,std::vector<std::string> addcomments);
#endif

  }; // end of class Grid_Creator_Base

} // end of namespace AMISIC

#include "Grid_Creator_Base.IC"
#include "Grid_Creator_Base.C"

#endif
