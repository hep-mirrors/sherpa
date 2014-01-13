#ifndef AMEGIC_Main_Single_Process_MHV_H
#define AMEGIC_Main_Single_Process_MHV_H

#include "AMEGIC++/Main/Process_Base.H"

#include "AMEGIC++/Main/Polarisation.H"
#include "AMEGIC++/Amplitude/Amplitude_Handler.H"
#include "ATOOLS/Math/Permutation.H"
// #include "ATOOLS/Phys/Spin_Correlation_Tensor.H"
#include <vector>
#include "AMEGIC++/Amplitude/Zfunctions/Mom.H"
#include "AMEGIC++/Amplitude/FullAmplitude_MHV_Base.H"

namespace AMEGIC {
  class Helicity;
  class Phase_Space_Generator;

  typedef std::set<std::pair<size_t,size_t> >     Combination_Set;
  typedef std::map<size_t,ATOOLS::Flavour_Vector> CFlavVector_Map;

  class Single_Process_MHV: public Process_Base {
  private:
    int                     m_gen_str;
    std::string             m_ptypename,m_libname,m_pslibname;
    Polarisation            m_pol;
    bool                    m_newlib;
    double                  m_iresult;
    int                     m_libnumb;
    bool                    m_ownamps;

    Helicity              * p_hel;
    Basic_Sfuncs          * p_BS;
    Amplitude_Handler     * p_ampl;
    String_Handler        * p_shand;
    Phase_Space_Generator * p_psgen;

    FullAmplitude_MHV_Base    * p_MHVamp;
    MomentumList              * p_momlist;
    Single_Process_MHV        * p_partner;

    Combination_Set m_ccombs;
    CFlavVector_Map m_cflavs;

    void FillCombinations(Point *const p,size_t &id);
    void FillCombinations();
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Process_MHV();
    ~Single_Process_MHV();
    
    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);

    bool Combinable(const size_t &idi,const size_t &idj);

    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);

    void           PolarizationNorm();
  private:

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  protected :
    int                 Tests();
    bool                CreateChannelLibrary();
  public:
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return m_newlib;}
    bool                FoundMappingFile(std::string &,std::string &);
    int                 InitAmplitude(MODEL::Model_Base *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    void                WriteMomFlavs(std::ofstream& to) { p_BS->WriteMomFlavs(to); }
    String_Handler    * GetStringHandler()             { return p_shand;}
    Amplitude_Handler * GetAmplitudeHandler()          { return p_ampl;}
    Helicity *          GetHelicity()                  { return p_hel; }    
    double              Result()                       { return m_iresult; } 

    int                             NumberOfDiagrams();
    AMEGIC::Point                 * Diagram(int i);
//     bool                            IsFreeOfFourVertex(AMEGIC::Point * _p);
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    std::string      LibName()                          { return m_libname;     }
    std::string      PSLibName()                        { return m_pslibname;   }
    AMEGIC::Process_Base   * Partner() const            { return p_partner;     }
    void             Minimize();

    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double         Partonic(const ATOOLS::Vec4D_Vector &,const int);
    double         DSigma(const ATOOLS::Vec4D_Vector &,bool);
    double         operator()(const ATOOLS::Vec4D *);
    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                        std::vector<std::vector<Complex> >& cols);
    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                        std::vector<std::vector<Complex> >& cols,
                        double sfactor);
    int PerformTests();

    int  Type() { return 1; }
  };
}



#endif
