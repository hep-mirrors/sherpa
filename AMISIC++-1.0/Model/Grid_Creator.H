#ifndef Grid_Creator_H
#define Grid_Creator_H

#include "XS_Group.H"
#include "Amisic_Histogram.H"
#include "File_IO_Base.H"

namespace AMISIC {

  class Simple_Chain;

  class Grid_Creator: public ATOOLS::File_IO_Base {
  public:
    
    typedef Amisic_Histogram<double> Amisic_Histogram_Type;
    typedef std::map<std::string,Amisic_Histogram_Type*> Amisic_Histogram_Map;
    
  private:
    
    Amisic_Histogram_Map *p_histograms;
    EXTRAXS::XS_Group    *p_processes;
    ATOOLS::Axis<double> *p_xaxis, *p_yaxis;

    const ATOOLS::Variable_Base<double> *p_variable;
    
    double m_gridxmin, m_gridxmax, m_griddeltax, m_binerror;

    std::string m_gridxvariable, m_gridxscaling;
    std::string m_gridyvariable, m_gridyscaling;
    
    std::string m_xsextension, m_mcextension, m_datatag;
    
    int m_criterion, m_outputlevel;

    long unsigned int m_initevents, m_maxevents, m_events;
    
    std::string MakeString(std::vector<std::string> input) const;
    
    bool CollectProcesses(EXTRAXS::XS_Base *const process);

    bool UpdateHistogram(EXTRAXS::XS_Base *const process);

    bool InitializeCalculation();
    bool CreateInitialGrid();
    bool CreateOptimizedGrid();
    
    friend class Simple_Chain;

  public:
    
    // constructor
    Grid_Creator(Amisic_Histogram_Map *histograms,
		 EXTRAXS::XS_Group *const processes);

    // destructor
    ~Grid_Creator();
    
    // member functions
    void Clear();

    bool ReadInGrid();
    bool CreateGrid();

    bool ReadInArguments(std::string tempifile=ATOOLS::nullstring,
			 std::string tempipath=ATOOLS::nullstring);
    bool WriteOutGrid(std::vector<std::string> addcomments=
		      std::vector<std::string>(0));

    // inline functions
    inline double GridXMin() const { return m_gridxmin; }
    inline double GridXMax() const { return m_gridxmax; }
   
    inline const std::string XSExtension() const  
    { return m_xsextension;  }
    inline const std::string MCExtension() const  
    { return m_mcextension;  }

    inline void SetXSExtension(const std::string xsextension)
    { m_xsextension=xsextension; }
    inline void SetMCExtension(const std::string mcextension)
    { m_mcextension=mcextension; }

    inline void SetGridXMin(const double gridxmin) 
    { m_gridxmin=gridxmin; }
    inline void SetGridXMax(const double gridxmax) 
    { m_gridxmax=gridxmax; }

  }; // end of class Grid_Creator

} // end of namespace AMISIC

#endif
