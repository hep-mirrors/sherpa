#ifndef AMEGIC_DipoleSubtraction_Single_OSTerm_H
#define AMEGIC_DipoleSubtraction_Single_OSTerm_H

#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "AMEGIC++/Main/Single_Process.H"

#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/Process_Info.H"
#include <vector>

namespace AMEGIC {
  class Helicity;
  class Amplitude_Handler;
  class DipoleSplitting_Base;
  class Single_OSTerm : public Process_Base {
  private:
    bool                    m_valid;
    Single_OSTerm     * p_partner;
    ATOOLS::Vec4D         * p_OS_mom;
    ATOOLS::Vec4D_Vector    p_OS_labmom;

    Process_Base           * p_os_process;
    PHASIC::Process_Info              m_osinfo;

    size_t                  m_pi,m_pj,m_pk;
    size_t                  m_switch;
//     size_t                  m_LOpij,m_LOpk;
    ATOOLS::Flavour         m_fli,m_flj,m_flij,m_flk;
    double                 m_wwindow;
    ATOOLS::NLO_subevt      m_subevt;
    
    PHASIC::Process_Integrator*     p_realint;

    std::vector<size_t> m_sids;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_OSTerm(const PHASIC::Process_Info &,size_t,size_t,size_t,PHASIC::Process_Integrator*);
    ~Single_OSTerm();

    void SetScale(const PHASIC::Scale_Setter_Arguments &args);

    /*------------------------------------------------------------------------------

      Generic stuff for initialization of Single_OSTermes

      ------------------------------------------------------------------------------*/

  private:
    void           BuildDecay(PHASIC::Subprocess_Info&);
    bool           DetermineType();
    std::vector<PHASIC::Subprocess_Info>::iterator
                   FindInInfo(PHASIC::Subprocess_Info& fi, int idx) const;
    void           SetLOMomenta(const ATOOLS::Vec4D*,const ATOOLS::Poincare &);
    double         lambda(double x, double y, double z);

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/

  public:
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return p_os_process->NewLibs();}
    int                 InitAmplitude(MODEL::Model_Base *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    String_Handler    * GetStringHandler()             {  return NULL;  }
    Amplitude_Handler * GetAmplitudeHandler()          { return p_partner->GetOSProcess()->GetAmplitudeHandler();}
    Helicity *          GetHelicity()                  {  return p_partner->GetOSProcess()->GetHelicity(); }    
    bool                IsValid()                      { return m_valid; }
    Process_Base*     GetOSProcess()                 { return p_os_process; }
    ATOOLS::NLO_subevt* GetSubevt()                    { return &m_subevt; }
    ATOOLS::Vec4D*      GetLOmom()                     { return p_OS_mom; }
    void                PrintProcessSummary(int it);

    int                 NumberOfDiagrams();
    AMEGIC::Point     * Diagram(int i);

    inline size_t       Li()                           { return m_pi; }
    inline size_t       Lj()                           { return m_pj; }
    inline size_t       Lk()                           { return m_pk; }
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    void             SetLookUp(const bool lookup);
    std::string      PSLibName()                        { return p_partner->GetOSProcess()->PSLibName();   }
    Process_Base   * Partner() const                    { return p_partner;     }
    void             Minimize();

    void SetSelector(const PHASIC::Selector_Key &key);
    
    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double         Partonic(const ATOOLS::Vec4D_Vector &,const int);
    double         operator()(const ATOOLS::Vec4D *,const ATOOLS::Poincare &cms,const int);

    void           FillAmplitudes(METOOLS::Amplitude_Tensor*, double =1.);

    int            Type() { return 15; }
    bool PerformTests() { return true; }
    void PrintLOmom();
  };
}



#endif

