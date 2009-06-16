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
    DipoleSplitting_Base*   p_dipole;
    PHASIC::NLO_subevt      m_subevt;
    double                  m_lastscale;
    
    PHASIC::Process_Integrator*     p_realint;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_DipoleTerm(const PHASIC::Process_Info &,size_t,size_t,size_t,PHASIC::Process_Integrator*);
    ~Single_DipoleTerm();

    /*------------------------------------------------------------------------------

      Generic stuff for initialization of Single_DipoleTermes

      ------------------------------------------------------------------------------*/
  private:
    bool           DetermineType();
    void           SetLOMomenta(const ATOOLS::Vec4D*);
    bool           CompareLOmom(const ATOOLS::Vec4D*);

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/

  public:
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return p_LO_process->NewLibs();}
    int                 InitAmplitude(MODEL::Model_Base *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    String_Handler    * GetStringHandler()             { return p_partner->GetLOProcess()->GetStringHandler();}
    Amplitude_Handler * GetAmplitudeHandler()          { return p_partner->GetLOProcess()->GetAmplitudeHandler();}
    Helicity *          GetHelicity()                  { return p_partner->GetLOProcess()->GetHelicity(); }    
    double              Result()                       { return m_iresult; } 
    bool                IsValid()                      { return m_valid; }
//     Vec4D*              GetLOMomenta()                 { return m_LO_mom; }
//     Flavour*            GetLOFlavs()                   { return m_LO_process->Flavour(); }
    Single_LOProcess*   GetLOProcess()                 { return p_LO_process; }
    dpt::dipoltype      GetDipoleType()                { return m_dipoletype; }
    PHASIC::NLO_subevt* GetSubevt()                    { return &m_subevt; }
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
    Process_Base   * Partner()                          { return p_partner;     }
    void             Minimize();
    
    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double         Differential(const ATOOLS::Vec4D_Vector &);
    double         Differential2();
    double         operator()(const ATOOLS::Vec4D *);

    double         GetDPSF(const ATOOLS::Vec4D *);

    void           FillAmplitudes(HELICITIES::Amplitude_Tensor*, double =1.);

    int            Type() { return 15; }
    bool PerformTests() { return true; }
    void PrintLOmom();
    const ATOOLS::Vec4D_Vector &ActiveMom() const { return p_LO_labmom; }
  };
}



#endif

