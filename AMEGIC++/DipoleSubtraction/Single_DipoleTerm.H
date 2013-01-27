#ifndef AMEGIC_DipoleSubtraction_Single_DipoleTerm_H
#define AMEGIC_DipoleSubtraction_Single_DipoleTerm_H

#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"

#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include <vector>

namespace AMEGIC {
  class Helicity;
  class Amplitude_Handler;
  class DipoleSplitting_Base;
  class Single_DipoleTerm : public Process_Base {
  private:
    double                  m_iresult;
    bool                    m_valid;
    double                  m_dalpha;
    double                  m_dkt2max;
    double                  m_maxgsmass;

    Single_DipoleTerm     * p_partner;

    Single_LOProcess      * p_LO_process;
    ATOOLS::Vec4D         * p_LO_mom;             
    ATOOLS::Vec4D_Vector    p_LO_labmom;

    dpt::dipoltype          m_dipoletype;
    int                     m_ftype;
    size_t                  m_pi,m_pj,m_pk;
    size_t                  m_LOpij,m_LOpk;
    ATOOLS::Flavour         m_fli,m_flj,m_flij,m_flk;
    ATOOLS::Flavour_Vector  m_lofl;
    DipoleSplitting_Base*   p_dipole;
    ATOOLS::NLO_subevt      m_subevt;
    double                  m_lastscale;
    
    PHASIC::Process_Integrator*     p_realint;

    std::vector<size_t> m_sids;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_DipoleTerm(const PHASIC::Process_Info &,size_t,size_t,size_t,PHASIC::Process_Integrator*);
    ~Single_DipoleTerm();

    void SetScale(const PHASIC::Scale_Setter_Arguments &args);
    void SetKFactor(const PHASIC::KFactor_Setter_Arguments &args);
    size_t SetMCMode(const size_t mcmode);
    size_t SetClusterMode(const size_t cmode);

    /*------------------------------------------------------------------------------

      Generic stuff for initialization of Single_DipoleTermes

      ------------------------------------------------------------------------------*/

  private:
    bool           DetermineType();
    std::vector<PHASIC::Subprocess_Info>::iterator
                   FindInInfo(PHASIC::Subprocess_Info& fi, int idx) const;
    void           SetLOMomenta(const ATOOLS::Vec4D*,const ATOOLS::Poincare &);
    bool           CompareLOmom(const ATOOLS::Vec4D*);

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/

  public:

    ATOOLS::Flavour ReMap(const ATOOLS::Flavour &fl,const size_t &id) const;
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return p_LO_process?p_LO_process->NewLibs():false;}
    int                 InitAmplitude(MODEL::Model_Base *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    bool                Trigger(const ATOOLS::Vec4D_Vector &p);
    String_Handler    * GetStringHandler()             { return p_partner->GetLOProcess()->GetStringHandler();}
    Amplitude_Handler * GetAmplitudeHandler()          { return p_partner->GetLOProcess()->GetAmplitudeHandler();}
    Helicity *          GetHelicity()                  { return p_partner->GetLOProcess()->GetHelicity(); }    
    double              Result()                       { return m_iresult; } 
    bool                IsValid()                      { return m_valid; }
//     Vec4D*              GetLOMomenta()                 { return m_LO_mom; }
//     Flavour*            GetLOFlavs()                   { return m_LO_process->Flavour(); }
    Single_LOProcess*   GetLOProcess()                 { return p_LO_process; }
    dpt::dipoltype      GetDipoleType()                { return m_dipoletype; }
    ATOOLS::NLO_subevt* GetSubevt()                    { return &m_subevt; }
    std::string         GetSplitConfID();
    ATOOLS::Vec4D*      GetLOmom()                     { return p_LO_mom; }
    void                PrintProcessSummary(int it);

    int                 NumberOfDiagrams();
    AMEGIC::Point     * Diagram(int i);

    inline size_t       Lijt()                         { return p_LO_process->GetEmit(); }
    inline size_t       Lkt()                          { return p_LO_process->GetSpect(); }
    inline size_t       Li()                           { return m_pi; }
    inline size_t       Lj()                           { return m_pj; }
    inline size_t       Lk()                           { return m_pk; }
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    void             SetLookUp(const bool lookup);
    std::string      PSLibName()                        { return p_partner->GetLOProcess()->PSLibName();   }
    Process_Base   * Partner() const                    { return p_partner;     }
    void             Minimize();

    void SetSelector(const PHASIC::Selector_Key &key);
    void SetShower(PDF::Shower_Base *const ps);
    
    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double         Partonic(const ATOOLS::Vec4D_Vector &,const int);
    double         operator()(const ATOOLS::Vec4D *,const ATOOLS::Poincare &cms,const int);

    int            Type() { return 15; }
    bool PerformTests() { return true; }
    void PrintLOmom();
  };
}



#endif

