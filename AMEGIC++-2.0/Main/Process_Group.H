#ifndef Process_Group_H
#define Process_Group_H

#include "Primitive_Observable_Base.H"
#include "Process_Base.H"
#include "Polarisation.H"
#include <string>

namespace AMEGIC {
  class Process_Group : public Process_Base {
  protected :
    std::vector<Process_Base *> m_procs;

    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public :
    Process_Group(int,int,APHYTOOLS::Flavour *&,
		  ISR::ISR_Handler *,BEAM::Beam_Spectra_Handler *,APHYTOOLS::Selector_Data *,
		  int,int,int,double,Pol_Info* _pl=0,int _nex=0,APHYTOOLS::Flavour * _ex_fl=0);
    Process_Group();
    ~Process_Group();
    /*------------------------------------------------------------------------------

      Management of Process_Groups

      ------------------------------------------------------------------------------*/
  private :
    void           ConstructProcesses(APHYTOOLS::Selector_Data *);
    void           GroupProcesses();
  public :
    void           Add(Process_Base *);
    void           SelectOne();
    void           DeSelect();
    Process_Base * Selected();
    void           Empty();
    int            Size()                          { return m_procs.size(); }
    Process_Base * operator[] (int idx)            { return m_procs[idx]; } 
    void           SetResDir(std::string);
    void           SetScale(double);
    void           SetISRThreshold(double);
    void           SetTables(bool);
    void           SetTotalXS(int);
    void           SetAtoms(bool);
    void           SetPartner(Single_Process * _p) { }
    Process_Base * Partner()                       { return 0; }
    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  public:
    int            InitAmplitude(Interaction_Model_Base *,Topology *,AMATOOLS::Vec4D *&,
				 std::vector<double> &,std::vector<Single_Process *> &);
    bool           SetUpIntegrator();
    void           InitAnalysis(std::vector<APHYTOOLS::Primitive_Observable_Base *>);
    /*------------------------------------------------------------------------------

      Calculating total cross sections and single event generation

      ------------------------------------------------------------------------------*/
  public:
    bool           CalculateTotalXSec();
    bool           LookUpXSec(double,bool,std::string);
    void           RescaleXSec(double);
    bool           PrepareXSecTables();
    void           AddPoint(const double value); 
    double         Differential(AMATOOLS::Vec4D *);
    double         Differential2();
    double         DSigma(AMATOOLS::Vec4D *,bool);
    double         DSigma2();
    bool           OneEvent();
    bool           SameEvent();
    double         WeightedEvent();
    /*------------------------------------------------------------------------------

      Helpers

      ------------------------------------------------------------------------------*/
  public :
    void           PrintDifferential();
    void           ControlOutput(AMATOOLS::Vec4D * moms);

  };
}

#endif


  /*!
    Organizes lists of processes as a list of Single_Processes and provides
    their integration channels as well. So, for integration we'll ALWAYS need
    the processes stuffed into Process_Groups. 

    There's one point related to that : We'll have to make sure, that only 
    identical numbers of F.S. particles with identical masses are organized in groups. 
    Otherwise we'll face desaster in the integration.
  */
    /*! 
      The list of single processes and the multichannel integrator of the group.
    */
    /*!
      In this constructor, some cuts and the Phase_Space_Hanlder are
      instantiated, all other parameters are set to start values.
    */
    /*!
      The empty constructor does basiscally nothing. It is just a wrapper to be filled
      later on. I wonder whether we really need it ...
    */
    /*!
      Add a process to the list of processes, provided that numbers of
      external legs etc. agree. This is to ensure that we can integrate
      a group in one go.
    */
    /*!
      InitAmplitude in the Process_Group iterates over all constituents and calls their
      respective InitAmplitude. If the group is atomic, i.e. if it is not to be integrated in
      one go, the momenta used for the checks are reset before each such call is made.
      This ensures that no unwanted cross-references are established.
      If all individually called InitAmplitudes yielded a "true" a "true" is returned.
    */ 
    /*!
      For Process_Groups, setting up an integrator is a little bit tricky.
      After the initialisation of the cuts and the multi channel for the final state
      part the latter will be set up by adding all channels of all processes that have
      no partner. With this condition the number of channels is at least a bit reduced to
      those that are potentially different. A similar thing happens for the 
      initial state part. Iterating over the same processes, types, masses and widths
      are added to corresponding vectors. After initializing the Phase_Space_Handler,
      the CreateIntegrator method will be called.
    */
    /*!
      This method steers the calculation of total cross sections.
      Depending on whether the group is atomar or not, either the methods 
      CalculateTotalXSEc of the processes are called or the following procedure
      is performed :
      First, cross references are established, see EstablishReferences().
      Finally the Phase_Space_Handler is employed to integrate the process via
      Integrate() and the result is stored internally via SetTotalXS():
    */

    /*!
      Adding the result at a point to the result(s) of the group and its processes.
      The value will be added to the totalsum of the group and, in due proportion,
      to the totalsum of every piece of the group.
    */
    /*!
      Stores the value of the total cross section in the internal variables and
      does the same for the single processes.
    */
