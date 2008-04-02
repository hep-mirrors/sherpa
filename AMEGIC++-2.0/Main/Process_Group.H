#ifndef Process_Group_H
#define Process_Group_H

#include "Process_Base.H"
#include "Polarisation.H"
#include <string>

#ifdef USING__Threading 
#include <pthread.h> 
#endif 

namespace AMEGIC {
  class Single_Process;

#ifdef USING__Threading
  class Process_Group;

  struct AME_PS_TID {
    pthread_t m_id;
    Process_Group *p_proc;
    double m_d;
    const ATOOLS::Vec4D *p_p;
    size_t m_s, m_m, m_b, m_e, m_i;
    pthread_mutex_t m_s_mtx, m_t_mtx;
    pthread_cond_t m_s_cnd, m_t_cnd;
    AME_PS_TID(Process_Group *const proc): 
      p_proc(proc), m_s(2), m_b(0), m_e(0) {}
  };// end of struct AME_PS_TID

  typedef std::vector<AME_PS_TID*> AME_PS_TID_Vector; 
#endif
 
  class Process_Group : public Process_Base {
  protected :
    std::vector<Process_Base *> m_procs;
    bool m_resetted;
    int  m_weventmode;
    int m_enable_mhv; 
    std::string m_mfname;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
#ifdef USING__Threading 
    std::vector<Process_Base *> m_umprocs, m_mprocs;
    AME_PS_TID_Vector m_cts;

    static void *TDSigma(void *arg);
#endif 
  public :
    Process_Group(Process_Info*,int,int,ATOOLS::Flavour *&,
		  PDF::ISR_Handler * =NULL,BEAM::Beam_Spectra_Handler * =NULL,
		  ATOOLS::Selector_Data * =NULL,
		  int=0,int=99,int=99,int=0,PHASIC::scl::scheme=PHASIC::scl::unknown,double=-1.,
		  Pol_Info* _pl=0,int _nex=0,ATOOLS::Flavour * _ex_fl=0, 
		  std::string cuttag="-1.", double error=-1.,
		  std::string e_func="1", int enable_mhv=0);
    Process_Group();
    ~Process_Group();
    /*------------------------------------------------------------------------------

      Management of Process_Groups

      ------------------------------------------------------------------------------*/
  private :
    int            SetPolarisations(char * plindex, Pol_Info * pl, int * beam_is_poled,const int &nout);
    void           ConstructProcesses(ATOOLS::Selector_Data *);
    void           PrepareTerminate();
    void           WriteMappingFile();
  public :
    bool           Find(std::string,Process_Base *&);
    void           Add(Process_Base *);
    void           GroupProcesses();
    bool           SelectOne();
    bool           SelectOneFromList();
    void           AddEvent(const double &xs);
    void           DeSelect();
    bool           ReSelect(int);
    void           AddEvent(const double xs,const double validxs,const int ncounts);
    void           SetEvents(const double number);
    void           ResetEvents(double gmin);
    void           GetGMin(double &g, double &meff);
    void           SetWEventMode(int mode);
    void           SetPSHandler(PHASIC::Phase_Space_Handler *const pshandler); 

    void           Empty();
    size_t         Size()                          { return m_procs.size(); }
    Process_Base * operator[] (int idx)            { return m_procs[idx]; } 
    void           SetResDir(std::string);
    void           SetMFname(std::string name)         { m_mfname=name; }
    void           SetISRThreshold(double);
    void           SetTables(bool);
    void           SetTotal(int flag, int depth=0);
    void           SetMax(const double max, int depth=0);
    void           ResetMax(int);
    void           SetMaxJetNumber(int max);
    void           SetCoreMaxJetNumber(int max);
    void           SetAtoms(bool);
    void           SetPartner(Single_Process * _p) { }
    Process_Base * Partner()                       { return 0; }
    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  public:
    int            InitAmplitude(MODEL::Model_Base *,Topology *,ATOOLS::Vec4D *&,
				 std::vector<Process_Base *> &,
				 std::vector<Process_Base *> &, int &, int &, int &);
    int            InitAmplitude(MODEL::Model_Base *,Topology *) { return 1; }
    bool           SetUpIntegrator();
    void           InitWeightHistogram();
    Amplitude_Handler * GetAmplitudeHandler() {return 0;} 
    Helicity *          GetHelicity() {return 0;}
    bool                NewLibs() {return 0;}
    std::string         PSLibName() {return std::string("");}        
    /*------------------------------------------------------------------------------

      Calculating total cross sections and single event generation

      ------------------------------------------------------------------------------*/
  public:
    bool           CalculateTotalXSec(std::string _resdir=std::string(""));
    bool           LookUpXSec(double,bool,std::string);
    void           RescaleXSec(double);
    void           SetupEnhance();
    void           SetFactorizationScale(const std::string &muf2);
    void           SetRenormalizationScale(const std::string &mur2);

    bool           PrepareXSecTables();
    void           WriteOutXSecs(std::ofstream &);
    void           WriteOutHistogram(std::string);
    void           AddPoint(const double value); 
    double         Differential(const ATOOLS::Vec4D *);
    double         Differential2();
    double         DSigma(const ATOOLS::Vec4D *,bool);
    double         DSigma2();
    ATOOLS::Blob_Data_Base * OneEvent(double = -1.);
    ATOOLS::Blob_Data_Base * SameEvent();
    ATOOLS::Blob_Data_Base * WeightedEventNS(const int mode=0);
    ATOOLS::Blob_Data_Base * WeightedEvent(const int mode=0);
    ATOOLS::Blob_Data_Base * SameWeightedEvent();
    void           OptimizeResult();
    /*------------------------------------------------------------------------------

      Helpers

      ------------------------------------------------------------------------------*/
  public :
    void           PrintDifferential();
    void           ControlOutput(ATOOLS::Vec4D * moms);
    void           SetPrintGraphs(bool print_graphs=true);

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
