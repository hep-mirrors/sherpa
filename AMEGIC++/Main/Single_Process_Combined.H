#ifndef AMEGIC_Main_Single_Process_Combined_H
#define AMEGIC_Main_Single_Process_Combined_H

#include "AMEGIC++/Main/Process_Base.H"

#include "AMEGIC++/Main/Polarisation.H"
#include "AMEGIC++/Amplitude/Amplitude_Handler.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
// #include "ATOOLS/Phys/Spin_Correlation_Tensor.H"
#include <vector>

namespace AMEGIC {
  class Helicity;
  class Phase_Space_Generator;

  typedef std::set<std::pair<size_t,size_t> >     Combination_Set;
  typedef std::map<size_t,ATOOLS::Flavour_Vector> CFlavVector_Map;

  class Single_Process_Combined: public Process_Base {
  private:
    int                     m_gen_str;
    std::string             m_ptypename,m_libname,m_pslibname;
    Polarisation            m_pol;
    bool                    m_newlib;
    double                  m_iresult;
    int                     m_libnumb;

    Helicity              * p_hel;
    Basic_Sfuncs          * p_BS;
    Amplitude_Handler     * p_ampl;
    String_Handler        * p_shand;
    Phase_Space_Generator * p_psgen;

    Single_Process_Combined * p_partner;

    PHASIC::Tree_ME2_Base * p_me2;

    Combination_Set m_ccombs;
    CFlavVector_Map m_cflavs;

    void FillCombinations(Point *const p,size_t &id);
    void FillCombinations();
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Process_Combined();
    ~Single_Process_Combined();
    
    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);

    bool Combinable(const size_t &idi,const size_t &idj);

    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);

  private:
    void           PolarizationNorm();

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  protected :
    int                 Tests();
    int                 TestLib();
    int                 CheckLibraries();
    int                 CheckStrings(Single_Process_Combined*);
    void                WriteLibrary();
    std::string         CreateLibName();
    void                CreateMappingFile(Single_Process_Combined*);
    void                UpdateMappingFile(std::string,std::map<std::string,Complex> &); 
    bool                CreateChannelLibrary();
    void                WriteAlternativeName(std::string aname);
    bool                CheckAlternatives(std::vector<Process_Base *> & links,
					  std::string procname);
  public:
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return m_newlib;}
    bool                FoundMappingFile(std::string &,std::string &);
    int                 InitAmplitude(MODEL::Model_Base *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    String_Handler    * GetStringHandler()             { return p_shand;}
    Amplitude_Handler * GetAmplitudeHandler()          { return p_ampl;}
    Helicity *          GetHelicity()                  { return p_hel; }    
    double              Result()                       { return m_iresult; } 

    int                             NumberOfDiagrams();
    AMEGIC::Point                 * Diagram(int i);
    bool                IsFreeOfFourVertex(Point * _p) { return 1-p_ampl->ExistFourVertex(_p); }
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    std::string      LibName()         { return m_libname;     }
    std::string      PSLibName()       { return m_pslibname;   }
    AMEGIC::Process_Base   * Partner() const { return p_partner;     }
    void             Minimize();

    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double         Partonic(const ATOOLS::Vec4D_Vector &,const int);
    double         DSigma(const ATOOLS::Vec4D_Vector &,bool);
    double         operator()(const ATOOLS::Vec4D *);
    int PerformTests();

  };
}



#endif
