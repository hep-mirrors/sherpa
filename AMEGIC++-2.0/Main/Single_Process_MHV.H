#ifndef Single_Process_MHV_H
#define Single_Process_MHV_H

#include "Process_Base.H"

#include "Polarisation.H"
#include "Helicity.H"
#include "Amplitude_Handler.H"
#include "Standard_Selector.H"
#include "../Amplitude/Zfunctions/MHV_Calculator.H"
#include "Permutation.H"
#include <vector>

namespace AMEGIC {
  class Single_Process_MHV : public Process_Base {
  private:
    std::string             m_ptypename,m_libname,m_pslibname;
    Polarisation            m_pol;
    double                  m_Norm;
    double                  m_save_max;
    bool                    m_newlib;
    double                  m_iresult;
    int                     m_libnumb;
    double                  m_sfactor;

    Helicity              * p_hel;
    Basic_Sfuncs          * p_BS;
    Amplitude_Handler     * p_ampl;
    String_Handler        * p_shand;
    MHV_Calculator        * p_MHVcalc;
    ATOOLS::Permutation   * p_perm;
    Single_Process_MHV    * p_partner;

//     double ** p_hcres;
//     double ** p_hcalpha;
//     int       m_scnt;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Process_MHV(int,int,ATOOLS::Flavour *,
		   PDF::ISR_Handler * =NULL,BEAM::Beam_Spectra_Handler * =NULL,
		   ATOOLS::Selector_Data * =NULL,
		   int=1,int=99,int=99,int=0,int=0,
		   double=-1.,Pol_Info* _pl=0,int _nex=0,ATOOLS::Flavour * _ex_fl=0,
		   int usepi=0,double ycut=-1.,double error=-1.,std::string e_func="1");
    Single_Process_MHV(Process_Info*,int,int,ATOOLS::Flavour *,
		   PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,
		   ATOOLS::Selector_Data *,
		   int,int,int,int,int,
		   double,Pol_Info* _pl,int _nex,ATOOLS::Flavour * _ex_fl,
		   int usepi,double ycut,double error,std::string e_func="1");
    ~Single_Process_MHV();
    void Add(Process_Base *) {
      ATOOLS::msg.Error()<<"Error in Single_Process_MHV. Method Add called."<<std::endl;
      abort();
    }

    /*------------------------------------------------------------------------------

      Generic stuff for initialization of Single_Process_MHVes

      ------------------------------------------------------------------------------*/
  private:
    void           PolarizationNorm();
    double         SymmetryFactors();
    void           PrepareTerminate();
    bool           CheckMapping(const Process_Base * proc);

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  protected :
    int                 Tests();
    int                 CheckLibraries();
    int                 CheckStrings(Single_Process_MHV*);
    void                WriteLibrary();
    std::string         CreateLibName();
    void                CreateMappingFile();
    bool                CreateChannelLibrary();
  public:
    bool                NewLibs() {return m_newlib;}
    bool                FoundMappingFile(std::string &,std::string &);
    int                 InitAmplitude(Interaction_Model_Base *,Topology *,ATOOLS::Vec4D *&,
				      std::vector<Process_Base *> &,std::vector<Process_Base *> &,
				      int &,int &,int &);
    int                 InitAmplitude(Interaction_Model_Base *,Topology *);
    void                InitDecay(Topology* top); 
    bool                SetUpIntegrator();
    bool                IsFreeOfFourVertex(Point * _p) { return 1-p_ampl->ExistFourVertex(_p); }
    String_Handler    * GetStringHandler()             { return p_shand;}
    Amplitude_Handler * GetAmplitudeHandler()          { return p_ampl;}
    Helicity *          GetHelicity()                  { return p_hel; }    
    double              Result()                       { return m_iresult; } 

    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    bool             SelectOne();
    bool             SelectOneFromList();
    void             DeSelect();
    bool             ReSelect(int);
    size_t           Size();
    bool             Find(std::string,Process_Base *&);
    void             WriteOutXSecs(std::ofstream & _to);
    Process_Base   * operator[] (int idx);
    std::string      LibName()                          { return m_libname;     }
    std::string      PSLibName()                        { return m_pslibname;   }
    void             SetTotal(int flag, int depth=0);
    void             FixISRThreshold();
    void             Empty();
    void             SetPartner(Single_Process_MHV * _p)    { p_partner = _p;       }
    Process_Base   * Partner()                          { return p_partner;     }
    void             Minimize();

    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    bool           CalculateTotalXSec(std::string _resdir=std::string(""));
    bool           LookUpXSec(double,bool,std::string);
    bool           PrepareXSecTables();
    void           AddPoint(const double value); 
    double         Differential(const ATOOLS::Vec4D *);
    double         Differential2();
    double         DSigma(const ATOOLS::Vec4D *,bool);
    double         DSigma2();
    double         operator()(const ATOOLS::Vec4D *);
    bool           OneEvent(double = -1.);
    bool           SameEvent();
    ATOOLS::Blob_Data_Base * WeightedEvent(const int mode=0);
    ATOOLS::Blob_Data_Base * WeightedEventNS(const int mode=0);
    ATOOLS::Blob_Data_Base * SameWeightedEvent();
    int            NumberOfDiagrams();
    Point        * Diagram(int i);
    void           ResetMax(int);
    void           OptimizeResult();

    /*------------------------------------------------------------------------------

      Helpers
      
      ------------------------------------------------------------------------------*/
  public:
    void           PrintDifferential();
    /*------------------------------------------------------------------------------

      Stuff to sample over helicities
      
      ------------------------------------------------------------------------------*/
  private:
    bool                    m_helsample,m_inithelsample;
    std::vector<int>        m_helnumbers;
    std::vector<double>     m_helalphas, m_helalphasaves;
    std::vector<double>     m_helresults, m_helresults2;
    std::vector<long int>   m_helthrows;
    long int                m_throws;
    int                     m_helnumber;
    double                  m_helresult, m_helresult2;

    void           InitializeHelicityWeights();
    const int      SelectedHelicity();
    void           AddToHelicity(const double,const int);
    void           OptimizeHelicityWeights();
    double         operator()(const ATOOLS::Vec4D *,const int);
  public:
    void           SetHelicitySampling(bool _helsample) { m_helsample = _helsample; }
  };
}



#endif

