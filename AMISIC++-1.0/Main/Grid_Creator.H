#ifndef Grid_Creator_H
#define Grid_Creator_H

#include "Grid_Handler.H"
#include "Simple_XS.H"
#include "Environment.H"
#include "Blob_List.H"
#include "Fixed_Variable_Channel.H"

namespace AMISIC {

  template <class Argument_Type,class Result_Type>
  class Grid_Creator: public Grid_Creator_Base<Argument_Type,Result_Type> {
  public:
    
    typedef Argument_Type GridArgumentType;
    typedef Result_Type GridResultType;
    
    typedef std::vector<GridArgumentType> GridArgumentVector;
    typedef std::vector<GridResultType> GridResultVector;

    typedef Grid_Handler<GridArgumentType,GridResultType> GridHandlerType;
    typedef std::vector<GridHandlerType*> GridHandlerVector;
    
    typedef Grid_Creator_Base<GridArgumentType,GridResultType> GridCreatorBaseType;
    typedef ATOOLS::Data_To_Function<GridArgumentType,GridResultType> GridFunctionType;
    
  private:
    
    GridHandlerVector   p_gridhandler;
    EXTRAXS::Simple_XS *p_processes;

    int  m_criterion;
    bool m_storemax;

    long int m_maxpoints[2];
    
    ATOOLS::Mom_Data m_initialdata;
    GridResultType m_lastxs, m_lastmax;
    
    ATOOLS::Axis<GridArgumentType> *p_xaxis;
    ATOOLS::Axis<GridResultType>   *p_yaxis;
    
    std::string m_xsextension, m_maxextension, m_mcextension;
    
  public:
    
    // constructor
    Grid_Creator(GridHandlerVector _p_gridhandler,EXTRAXS::Simple_XS *_p_processes);

    // destructor
    ~Grid_Creator();
    
    // setters
    void SetXSExtension(const std::string _m_xsextension);
    void SetMaxExtension(const std::string _m_maxextension);
    void SetMCExtension(const std::string _m_mcextension);

    void SetMaxPoints(const long int _m_maxpoints);
    
    // getters
    const std::string XSExtension() const;
    const std::string MaxExtension() const;
    const std::string MCExtension() const;

    const long int MaxPoints() const;
    
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
    
  }; // end of class Simple_Chain::Grid_Creator

} // end of namespace AMISIC

#include "Grid_Creator.inl.H"
#include "Grid_Creator.C"

#endif
