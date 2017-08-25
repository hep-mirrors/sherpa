#ifndef AMEGIC_DipoleSubtraction_Single_Virtual_Correction_H
#define AMEGIC_DipoleSubtraction_Single_Virtual_Correction_H

#include "AMEGIC++/Main/Process_Base.H"

#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"
#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "PHASIC++/Process/KP_Terms.H"

namespace PDF { class PDF_Base; } 

namespace PHASIC {
  class Virtual_ME2_Base;
  class Multi_Channel;
}

namespace AMEGIC {
  class Helicity;
  class Amplitude_Handler;
  class Phase_Space_Generator;

  class Single_Virtual_Correction : public Process_Base {
  protected:
    std::string             m_ptypename,m_pslibname;
    ATOOLS::sbt::subtype    m_stype;
    ATOOLS::sbt::subtype    m_ipart;
    ATOOLS::cs_itype::type  m_itype;
    double                  m_iresult;
    Phase_Space_Generator * p_psgen;
    PHASIC::Multi_Channel * p_fsmc;

    std::vector<std::vector<double> > m_dsijqcd,m_dsijew;
    std::vector<std::vector<double> > m_Q2ij;

    Single_Virtual_Correction        * p_partner;
    Single_LOProcess                 * p_LO_process;
    PHASIC::Massive_Kernels          * p_kernel_qcd, * p_kernel_ew;
    PHASIC::KP_Terms                 * p_kpterms_qcd, * p_kpterms_ew;
    PHASIC::Virtual_ME2_Base         * p_loopme;

    double *                p_reqborn;
    double                  m_cmur[2];
    double                  m_x0,m_x1,m_eta0,m_eta1,m_z0,m_z1;

    bool   m_force_init, m_sccmur, m_murcoeffvirt, m_loopmapped;
    bool   m_checkborn, m_checkpoles, m_checkfinite;
    size_t m_pspisrecscheme, m_pspfsrecscheme;
    size_t m_pspissplscheme, m_pspfssplscheme;
    double m_checkthreshold;
    int    m_bvimode, m_user_bvimode, m_epsmode, m_drmode, m_checkloopmap;

    double m_bsum, m_vsum, m_isum, m_n;
    double m_mbsum, m_mvsum, m_misum, m_mn;
    double m_lastb, m_lastv, m_lasti, m_lastkp, m_lastki;
    double m_finite, m_singlepole, m_doublepole;

    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Virtual_Correction();
    ~Single_Virtual_Correction();

    void AddPoint(const double &value);
    bool ReadIn(const std::string &pid);
    void WriteOut(const std::string &pid);
    void MPICollect(std::vector<double> &sv,size_t &i);
    void MPIReturn(std::vector<double> &sv,size_t &i);
    void MPISync(const int mode=0);

    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);
    void SetScale(const PHASIC::Scale_Setter_Arguments &args);
    void SetVariationWeights(ATOOLS::Variation_Weights *const vw);

    void SetSelector(const PHASIC::Selector_Key &key);
    void SetShower(PDF::Shower_Base *const ps);
    void SetNLOMC(PDF::NLOMC_Base *const mc);
    void SetFixedScale(const std::vector<double> &s);
    void SetSelectorOn(const bool on);
    void SetGenerator(PHASIC::ME_Generator_Base *const gen);
    ATOOLS::Flavour ReMap(const ATOOLS::Flavour &fl,const size_t &id) const;
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
    int                 InitAmplitude(Amegic_Model *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    void                RequestVariables(PHASIC::Phase_Space_Handler *const psh);
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
    double         operator()(const ATOOLS::Vec4D_Vector &,const int mode);
    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                        std::vector<std::vector<Complex> >& cols);

    int            NumberOfDiagrams();
    Point        * Diagram(int i);

    bool Combinable(const size_t &idi,const size_t &idj);
    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);
    void FillProcessMap(PHASIC::NLOTypeStringProcessMap_Map *apmap);

    inline double * RequestedBorn() { return p_reqborn; }

    inline const PHASIC::KP_Terms*         KPTermsQCD() { return p_kpterms_qcd; }
    inline const PHASIC::KP_Terms*         KPTermsEW()  { return p_kpterms_ew; }
    inline const PHASIC::Massive_Kernels * KernelQCD()  { return p_kernel_qcd; }
    inline const PHASIC::Massive_Kernels * KernelEW()   { return p_kernel_ew; }

    inline double Finite()     const { return m_finite; }
    inline double SinglePole() const { return m_singlepole; }
    inline double DoublePole() const { return m_doublepole; }

    inline void SetKPz0(const double &z0) { m_z0=z0; }
    inline void SetKPz1(const double &z1) { m_z1=z1; }

  protected:
    double         Calc_B();
    double         Calc_V(const ATOOLS::Vec4D_Vector&);
    double         Calc_V_WhenMapped(const ATOOLS::Vec4D_Vector&);
    double         Calc_I(const ATOOLS::Vec4D_Vector&);
    double         Calc_I(const ATOOLS::sbt::subtype,
                          const std::vector<size_t>&,
                          PHASIC::Massive_Kernels *,
                          const PHASIC::KP_Terms*,
                          const ATOOLS::Vec4D_Vector &,
                          std::vector<std::vector<double> >&);
    void           Calc_KP(const ATOOLS::Vec4D_Vector&);
    double         KPTerms(int mode, PDF::PDF_Base *pdfa,
                                     PDF::PDF_Base *pdfb, double scalefac2=1.0);
    double         Get_KPTerms(PDF::PDF_Base *pdfa, PDF::PDF_Base *pdfb,
                               const double&,const double&,
                               const ATOOLS::Flavour&,const ATOOLS::Flavour&,
                               const double&);
    void           CheckBorn();
    void           CheckFinite(const double &, const double &);
    void           CheckPoleCancelation(const ATOOLS::Vec4D_Vector);
    void           FillMEwgts(ATOOLS::ME_Weight_Info&);
    void           EndOptimize();
    void           AttachChargeFactors();
    void           ComputeChargeFactors();
    ATOOLS::ist::itype  AssignType(const size_t& id,
                                   const ATOOLS::sbt::subtype);
    bool           AllowAsSpecInPFF(const size_t& i,const size_t& k);

    void           PrintDSij();
  };
}



#endif

