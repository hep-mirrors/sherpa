#ifndef AMISIC_Model_Grid_Creator_H
#define AMISIC_Model_Grid_Creator_H
/*!
  \file Grid_Creator.H
  \brief Declares the class Grid_Creator
*/

#include "EXTRA_XS/Main/Process_Group.H"
#include "AMISIC++/Tools/Amisic_Histogram.H"
#include "ATOOLS/Org/File_IO_Base.H"

namespace AMISIC {

  class Simple_Chain;

  class Grid_Creator: public ATOOLS::File_IO_Base {
  public:
    
    typedef Amisic_Histogram<double> Amisic_Histogram_Type;
    typedef std::map<std::string,Amisic_Histogram_Type*> Amisic_Histogram_Map;
    
  private:
    
    Amisic_Histogram_Map *p_histograms;
    EXTRAXS::Process_Group    *p_processes;
    ATOOLS::Axis<double> *p_xaxis, *p_yaxis;

    const ATOOLS::Variable_Base<double> *p_variable;
    
    double m_gridxmin, m_gridxmax, m_griddeltax, m_binerror;

    std::string m_gridxvariable, m_gridxscaling;
    std::string m_gridyvariable, m_gridyscaling;
    
    std::string m_xsextension, m_mcextension, m_datatag, m_criterion;
    
    int m_outputlevel;

    long unsigned int m_initevents, m_maxevents, m_events;
    
    std::string MakeString(std::vector<std::string> input) const;
    
    bool CollectProcesses(PHASIC::Process_Base *const process);
    /*!
      \fn bool CollectProcesses(EXTRAXS::XS_Base *const process)
      \brief Examines the given XS_Base and extracts Single_Process'.
      Creates a histogram for each Single_Process and stores the histogram 
      pointer to the histogram map under the respective process name.
    */

    bool UpdateHistogram(PHASIC::Process_Base *const process);
    /*!
      \fn bool UpdateHistogram(EXTRAXS::XS_Base *const process)
      \brief Iterates over the Single_Process' contained in the given 
      process, determines their differential cross sections and 
      fills the respective histograms.
    */

    bool CreateGrid(PHASIC::Process_Group *const process);
    /*!
      \fn bool CreateInitialGrid()
      \brief Creates the grid.
    */
    
    friend class Simple_Chain;

  public:
    
    // constructor
    Grid_Creator(Amisic_Histogram_Map *histograms,
		 EXTRAXS::Process_Group *const processes);
    /*!
      \fn Grid_Creator(Amisic_Histogram_Map *histograms,
      EXTRAXS::Process_Group *const processes)
      \brief The standard constructor. For each of the processes
      contained in the given group Grid_Creator creates a histogram.
      The respective pointer is stored in the given map 
      under the name of the respective process.
    */

    // destructor
    ~Grid_Creator();
    
    // member functions
    void Clear();
    /*!
      \fn void Clear()
      \brief Clears the Grid_Creator.
    */

    bool ReadInGrid();
    /*!
      \fn bool ReadInGrid()
      \brief Reads the grids from exernal files. If the reading fails
      or the read grids have too less points, the method returns false.
    */
    bool CreateGrid();
    /*!
      \fn bool CreateGrid()
      \brief Creates the grids. 
    */

    bool InitializeCalculation(PHASIC::Process_Group *const process);
    /*!
      \fn bool InitializeCalculation()
      \brief Initializes the calculation.
    */

    bool ReadInArguments(std::string tempifile=ATOOLS::nullstring,
			 std::string tempipath=ATOOLS::nullstring);
    /*!
      \fn bool ReadInArguments(std::string tempifile=ATOOLS::nullstring,
      std::string tempipath=ATOOLS::nullstring)
      \brief Reads the basic information from the given file.
    */
    bool WriteOutGrid(std::vector<std::string> addcomments=
		      std::vector<std::string>(0),
		      const std::string &path="");
    /*!
      \fn bool WriteOutGrid(std::vector<std::string> addcomments=
      std::vector<std::string>(0),const std::string &path)
      \brief Writes the created grids to external files which are named
      according to the respective process.
    */

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
  /*!
    \class Grid_Creator
    \brief Creates histograms of variable spectra
    
    This class is employed to create histograms of variable 
    spectra, with special emphasis on transverse momenta. 
    It receives a list of processes from the internal process 
    library which is initialized and used to generate weighted events.
    During the event generation the corresponding histograms are filled.
  */

} // end of namespace AMISIC

#endif
