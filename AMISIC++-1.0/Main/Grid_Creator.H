#ifndef Grid_Creator_H
#define Grid_Creator_H

#include "Grid_Handler.H"
#include "XS_Group.H"
#include "Amisic_Histogram.H"

namespace AMISIC {

  class Simple_Chain;

  class Grid_Creator: public ATOOLS::File_IO_Base {
  public:
    
    typedef Amisic_Histogram<double>                Amisic_Histogram_Type;
    typedef Grid_Handler<double,double>             Grid_Handler_Type;
    typedef ATOOLS::Data_To_Function<double,double> Grid_Function_Type;

    typedef std::map<std::string,Grid_Handler_Type*>     Grid_Handler_Map;
    typedef std::map<std::string,Amisic_Histogram_Type*> Amisic_Histogram_Map;

    typedef std::vector<Grid_Handler_Type*> Grid_Handler_Vector;
    
  private:
    
    Grid_Handler_Map     *p_gridhandlers, *p_maxhandlers;
    Amisic_Histogram_Map  m_histograms;

    EXTRAXS::XS_Group *p_processes;
    ATOOLS::Mom_Data   m_momdata;
   
    ATOOLS::Axis<double> *p_xaxis, *p_yaxis;
    
    double m_gridxmin, m_gridxmax, m_griddeltax, m_binerror;

    std::string m_gridxvariable, m_gridxscaling;
    std::string m_gridyvariable, m_gridyscaling;
    
    std::string m_xsextension, m_maxextension, m_mcextension;
    
    bool m_storemax;
    int  m_criterion, m_outputlevel;

    long unsigned int m_initevents, m_maxevents, m_events;
    
    const ATOOLS::Vec4D *p_momenta;
    
    std::string MakeString(std::vector<std::string> input) const;
    
    bool CollectProcesses(EXTRAXS::XS_Base *const process);

    bool UpdateHistogram(EXTRAXS::XS_Base *const process,
			 const size_t trials);

    bool ExportHistogram(const std::string &name) const;

    bool InitializeCalculation();
    bool CreateInitialGrid();
    bool CreateOptimizedGrid();
    bool TranslateGrid();
    
    bool WriteSingleGrid(Grid_Handler_Type *grid,
			 std::vector<std::string> addcomments);

    friend class Simple_Chain;

  public:
    
    // constructor
    Grid_Creator(Grid_Handler_Map *gridhandlers,
		 Grid_Handler_Map *maxhandlers,
		 EXTRAXS::XS_Group *const processes);

    // destructor
    ~Grid_Creator();
    
    // member functions
    bool ReadInGrid();
    bool CreateGrid();

    bool ReadInArguments(std::string tempifile=ATOOLS::nullstring,
			 std::string tempipath=ATOOLS::nullstring);
    bool WriteOutGrid(std::vector<std::string> addcomments=
		      std::vector<std::string>(0),
		      std::string tempopath=ATOOLS::nullstring);

    // inline functions
    inline const std::string XSExtension() const  { return m_xsextension;  }
    inline const std::string MaxExtension() const { return m_maxextension; }
    inline const std::string MCExtension() const  { return m_mcextension;  }

    inline void SetXSExtension(const std::string xsextension)
    { m_xsextension=xsextension; }
    inline void SetMaxExtension(const std::string maxextension)
    { m_maxextension=maxextension; }
    inline void SetMCExtension(const std::string mcextension)
    { m_mcextension=mcextension; }

  }; // end of class Grid_Creator

} // end of namespace AMISIC

#endif
