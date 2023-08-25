#ifndef YFS_Main_YFS_Handler_H
#define YFS_Main_YFS_Handler_H

#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/ISR.H"
#include "YFS/Main/FSR.H"
#include "YFS/Main/Coulomb.H"
#include "YFS/Main/Semi.H"
#include "YFS/Tools/Debug.H"
#include "YFS/NLO/Virtual_ff.H"
#include "YFS/NLO/EEX.H"
#include "YFS/NLO/NLO_Base.H"
#include "BEAM/Main/Beam_Base.H"
#include "YFS/Main/Semi.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/CEEX/Ceex_Base.H"
#include "YFS/NLO/NLO_Base.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"



#include <iostream>
#include <fstream>
namespace ATOOLS { class Blob_Data_Base; }
// namespace BEAM   { class Beam_Base;      }
namespace YFS    {class Real_ff;         }
namespace YFS    {class Coulomb;         }
// namespace YFS    {class ISR;         }
// namespace PHASIC {class Tree_ME2_Base;         }
namespace OpenLoops {class OpenLoops_Interface; }


namespace YFS {
class Semi_YFS;
typedef std::vector<double> Double_Vector;


// class Remnant_Base;
class YFS_Handler : public YFS_Base {
    YFS::Virtual_ff      *p_virtff;
    YFS::Semi_YFS *SY;


public:
     // constructor
    YFS_Handler();

    // destructor
    ~YFS_Handler();
    YFS::Real_ff         *p_realff;
    // double m_mass2[2], m_x[1], m_u;
    double m_splimits[3];
    // double m_fixed_smin, m_fixed_smax, m_s;
    double m_spp, m_tchannelcut;
    double m_yfsweight, m_form, m_beta1, m_fsrW, m_isrPhaseSpaceW, m_fsrPhaseSpaceW;
    bool   m_bornInit, m_born_set, m_setparticles, m_asymbeams;
    double m_real, m_virt, m_yfsME, m_corrwgt, m_cms, m_oneloop;
    double m_nbar, m_NbarFSR, m_the, m_phi, m_cth, m_formfactor, m_betatWW, m_coulSub, m_ww_formfact;
    // int m_yfsmode, m_N, m_order, m_nmax, m_betaorder, m_N_FSR, m_fullform, m_formWW;
    // int m_fsrmode, m_useint, m_check_mass_reg, m_check_poles, m_use_fsr_beta;
    YFS::Coulomb *p_coulomb;
    PHASIC::Process_Base *proc;
    ATOOLS::Particle_Vector m_particles,  dipoleISR, dipoleFSR;
    ATOOLS::ParticleMomMap m_inparticles, m_outparticles;

    // std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR;
    ATOOLS::Poincare m_cmsboost;
    ATOOLS::Vec4D    m_fixvecs[2], p_cms[2], m_realPho, m_Q, m_P, m_PX, m_Wm, m_Wp;
    ATOOLS::Flavour m_bunch;
    std::vector<double> m_info_lab, m_info_cms, m_dvir, m_y, m_z, m_eikonalsISR, m_eikonalsFSR, m_realcorr;
    // Remnant_Base    *p_remnants[2];
    Vec4D_Vector m_plab, m_reallab, m_FSRPhotons, m_ISRPhotons, m_bornMomenta, m_momeikonal;
    ATOOLS::Flavour_Vector m_flavs;
    std::vector<double> m_mass;
    ATOOLS::Momenta_Stretcher m_stretcher;


    double m_f;


    bool CheckRemnantKinematics(const ATOOLS::Flavour &flav,
                                double &x, int beam, bool swaped);

public:

    YFS::ISR *p_isr;
    YFS::Debug *p_debug;
    // int m_mode, m_rmode, m_semiyfs, m_NReal, m_qedmode;
    ATOOLS::Vec4D_Vector m_p, m_beams;
    ATOOLS::nlo_type::code m_nlotype;
    std::vector<double> m_MEorder;
    BEAM::Beam_Base *p_beams[2];
    // YFS::Weight *p_yfsweight;
    YFS::Define_Dipoles *p_dipole;
    Vec4D m_photonSumISR, m_photonSumFSR, m_kboostVec, m_K, m_QN, m_beamP, m_beamM;
    Vec4D_Vector m_finalFSRPhotons, m_FSRMomentaME, m_fsrphotonsforME;
    Flavour_Vector m_inflav;
    Flavour_Vector m_outflav;
    EXTRAXS::ME2_Base *p_me;
    PHASIC::Process_Info p_pi;
    YFS::NLO_Base *p_nlo;
    // YFS::Virtual  *p_virt;
    // YFS::Real  *p_real;
    // YFS::RealReal  *p_realreal;
    // YFS::RealVirtual *p_realvirt;
    MODEL::Coupling_Map  m_cpls;
    MODEL::Model_Base  *p_model;
    YFS::FSR *p_fsr;
    YFS::YFS_Form_Factor *p_yfsFormFact;
    YFS::Define_Dipoles *p_dipoles;
    YFS::Ceex_Base *p_ceex;
    YFS::Dipole_Vector m_olddipoles, m_newdipoles;
    Poincare *p_boost, *p_Qboost, *p_Pboost, *p_rotate;
   

    // member functions
    void Init(std::vector<double> &splimits, const Flavour_Vector &fl);

    void SetYFSMember() const;
    // Structure_Function(PDF::PDF_Base *,ATOOLS::Flavour) *p_sf;
//
    bool On();
    bool MakeYFS(ATOOLS::Vec4D_Vector &p);
    void MakeYFS();
    void SetSprimeLimits(std::vector<double> &splimits);
    void SetFlavours(const Flavour_Vector &fl);
    void SetMomenta(const ATOOLS::Vec4D_Vector &plab);
    void SetBornMomenta(const ATOOLS::Vec4D_Vector &plab);
    void InitializeCEEX(const ATOOLS::Flavour_Vector &fl);
    void CalculateBeta();
    void GenerateWeight();
    void CreatMomentumMap();
    void MakeCEEX(const Vec4D_Vector &p);
    Vec4D_Vector GetPhotons();
    double Flux(const Vec4D&, const Vec4D&);
    inline double GetWeight() const      { return m_yfsweight;  }
    inline ATOOLS::Vec4D GetPhotonSum() const      { return m_photonSumISR;  }
    inline ATOOLS::Vec4D GetPhotonSumFSR() const      { return m_photonSumFSR;  }
    inline ATOOLS::Vec4D_Vector GetISRPhotons() const      { return m_ISRPhotons;  }
    inline ATOOLS::Vec4D_Vector GetFSRPhotons() const      { return m_FSRPhotons;  }
    // inline YFS::Real_ff GetNLO() {return p_realff;}
    inline void SetBorn(double b) { m_born = b; }
    inline double GetBorn() { return m_born; }
    inline double GetBeta1() const { return m_beta1; }
    inline double GetReal() const { return m_real; }
    inline double GetVirtual() const { return m_virt; }
    inline double GetV()   const { return m_v;}
    inline void SetMode(int mode)    { m_mode = mode;}
    inline void SetModel(MODEL::Model_Base *model)    { p_model = model;}
    inline void SetOrder(int order)  { m_order = order;}
    inline void SetProcessInfo(PHASIC::Process_Info &pi) {p_pi = pi;}
    inline void SetISRPhotons(Vec4D_Vector &k) {m_ISRPhotons = k;}
    inline int  GetMode()   const {return m_mode; }
    inline int  GetBetaOrder()   const {return m_betaorder; }
    inline int  GetOrder()   const {return m_order; }
    inline int  GetFSRMode()   const {return m_fsrmode; }
    inline int  NPhoton()   const {return m_N; }
    inline bool UsingLoopTool()   const {return m_looptool; }
    inline bool FillBlob()   const {return m_fillblob; }
    inline bool CalculatingReal()   const {return m_YFS_RS; }
    inline Vec4D_Vector GetMomenta() { return m_plab; }
    inline Vec4D GetMomenta(int i) { return m_plab[i]; }
    inline Vec4D_Vector BornMomenta() { return m_bornMomenta; }
    inline int  GetRealPhotonN()    const {return m_NReal;}
    inline void SetRealPhotonN(int n)  {m_NReal = n;}
    inline void SetQedMode(int qed)  {m_qedmode = qed;}
    inline int  GetQedMode()  {return m_qedmode;}
    inline void SetSprimeMin(double sp) {m_smin = sp;}
    inline void SetCouplings(MODEL::Coupling_Map &cpls)   { m_cpls = cpls; }
    inline void SetNLOType(ATOOLS::nlo_type::code nlo_type) {m_nlotype = nlo_type;}
    inline NLO_Base *NLO() {return p_nlo; }
    // inline void SetFSRMomenta(ATOOLS::Vec4D_Vector f) { m_FSRMomentaME = f;}
    inline std::vector<double> SetOrder(double qcd, double ew) {
        m_MEorder.clear();
        m_MEorder.push_back(qcd);
        m_MEorder.push_back(ew);
        return m_MEorder;
    } ;
    void SetBeam(BEAM::Beam_Spectra_Handler *beam);
    // void InitializeVirtual(const PHASIC::Process_Info& pi);
    // void InitializeReal(const PHASIC::Process_Info& pi);
    // void InitializeRealReal(const PHASIC::Process_Info& pi);
    // void InitializeRealVirtual(const PHASIC::Process_Info& pi);
    // double CalculateVirtual(Vec4D_Vector p, double);
    // double CalculateReal();
    // double CalculateRealVirtual();
    // double CalculateRealReal();

    inline void SetSprime(double sp) { m_sp = sp; }
    inline double GetPhaseSpaceWeight() { return m_isrPhaseSpaceW; }
    inline void SetPhaseSpaceWeight(double wgt) { m_fsrPhaseSpaceW = wgt; }
    inline void SetRunMode(const int &rmode) { m_rmode = rmode; }
    inline bool UseCEEX() {return m_useceex;}
    inline YFS::Ceex_Base *CEEX() {return p_ceex;}

    double Eikonal(Vec4D k);
    double CalculateNLO();
    double MasslessEikonal(Vec4D k, Vec4D p1, Vec4D p2);
    bool CalculateFSR();
    bool CalculateISR();
    bool CalculateFSR(Vec4D_Vector &k);
    bool CheckMomentumConservation();
    void CheckMasses();
    void Reset();
    int NHardPhotons(const Vec4D_Vector &k);
    double Flux();
    void AddFormFactor();
    void MakeWWVecs(ATOOLS::Vec4D_Vector p);
    void CalculateCoulomb();
    void CalculateWWForm();
};// end of class YFS_Handler



}// end of namespace YFS

#endif
