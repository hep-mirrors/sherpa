#ifndef AMEGIC_DipoleSubtraction_Single_Virtual_Correction_H
#define AMEGIC_DipoleSubtraction_Single_Virtual_Correction_H

#include "AMEGIC++/Main/Process_Base.H"

#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include <vector>

namespace PDF { class PDF_Base; } 

namespace PHASIC { class Virtual_ME2_Base; }

namespace AMEGIC {
  class Helicity;
  class Amplitude_Handler;
  class Phase_Space_Generator;
  class DipoleSplitting_Base;
  class Flavour_KernelsA;
  class Massive_KernelsA;

  struct xp_coeff {
    double xp;
    double kpc;
  };

  class Single_Virtual_Correction : public Process_Base {
  protected:
    std::string             m_ptypename,m_pslibname;
    double                  m_iresult;
    double                  m_dalpha,m_dalpha_ff,m_dalpha_fi,m_dalpha_if,m_dalpha_ii;
    Phase_Space_Generator * p_psgen;

    double               ** p_dsij;
    double                  m_kpca[8];
    double                  m_kpcb[8];
    double                  m_cmur[2];
    double                  m_x0,m_x1;
    std::vector<xp_coeff>   m_xpa, m_xpb;

    Single_Virtual_Correction        * p_partner;

    Single_LOProcess      * p_LO_process;
    DipoleSplitting_Base  * p_dipole;
    Flavour_KernelsA      * p_flkern;
    Massive_KernelsA      * p_masskern;
    PHASIC::Virtual_ME2_Base* p_loopme;  

    bool m_massive, m_force_init;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Virtual_Correction();
    ~Single_Virtual_Correction();

    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);
    void SetScale(const PHASIC::Scale_Setter_Arguments &args);

    void SetSelector(const PHASIC::Selector_Key &key);
    void SetShower(PDF::Shower_Base *const ps);
    void SetFixedScale(const std::vector<double> &s);
    void SetSelectorOn(const bool on);
    void SetGenerator(PHASIC::ME_Generator_Base *const gen);
    /*------------------------------------------------------------------------------

      Generic stuff for initialization of Single_Virtual_Correctiones

      ------------------------------------------------------------------------------*/
  protected:
    void                PolarizationNorm();
    double              Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    bool                CreateChannelLibrary();

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  public:
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs();
    int                 InitAmplitude(MODEL::Model_Base *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    Single_LOProcess*   GetLOProcess()                 { return p_LO_process; }
    Amplitude_Handler * GetAmplitudeHandler()          { return p_partner->GetLOProcess()->GetAmplitudeHandler();}
    Helicity *          GetHelicity()                  { return p_partner->GetLOProcess()->GetHelicity(); }    
    double              Result()                       { return m_iresult; } 

    void                SelectLoopProcess();
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    void             SetLookUp(const bool lookup);
    std::string      LibName()                          { return p_partner->GetLOProcess()->LibName();     }
    std::string      PSLibName()                        { return p_partner->GetLOProcess()->PSLibName();   }
    Process_Base   * Partner()       const              { return p_partner;     }
    void             Minimize();

    int Type() { return 200; }
    
    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double         Partonic(const ATOOLS::Vec4D_Vector &,const int mode);
    double         DSigma(const ATOOLS::Vec4D_Vector &,bool,const int);
    double         DSigma2();
    double         operator()(const ATOOLS::Vec4D_Vector &,const int mode);

    void           FillAmplitudes(METOOLS::Amplitude_Tensor*, double =1.);
    int            NumberOfDiagrams();
    Point        * Diagram(int i);

    bool PerformTests() { return true; }

    bool Combinable(const size_t &idi,const size_t &idj);
    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);

  protected:
    double         Calc_Imassive(const ATOOLS::Vec4D *mom);
    double         Calc_I(const ATOOLS::Vec4D*);
    void           Calc_KP(const ATOOLS::Vec4D*, double, double, double, double, double);
    double         Get_KPterms(PDF::PDF_Base *pdfa, PDF::PDF_Base *pdfb,
			       const double&,const double&,ATOOLS::Flavour_Vector&); 
    void           CheckPoleCancelation(const ATOOLS::Vec4D*);
    void           FillMEwgts(ATOOLS::ME_wgtinfo&);

  };
}



#endif

