#ifndef YFS_Main_YFS_Handler_H
#define YFS_Main_YFS_Handler_H

#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/ISR.H"
#include "YFS/Main/FSR.H"
#include "YFS/Main/Coulomb.H"
#include "YFS/Tools/Debug.H"
#include "YFS/NLO/NLO_Base.H"
#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/CEEX/Ceex_Base.H"
#include "YFS/NLO/NLO_Base.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"



#include <iostream>
#include <fstream>
namespace ATOOLS { class Blob_Data_Base; }
namespace YFS    {class Real_ff;         }
namespace YFS    {class Coulomb;         }


namespace YFS {
typedef std::vector<double> Double_Vector;

class YFS_Handler : public YFS_Base {

public:
     // constructor
    YFS_Handler();

    // destructor
    ~YFS_Handler();
    double m_splimits[3];
    double m_yfsweight, m_form, m_fsrW;
    bool   m_setparticles, m_asymbeams;
    double m_real, m_virt;
    double m_nbar, m_formfactor, m_betatWW, m_coulSub, m_ww_formfact;

    YFS::Coulomb *p_coulomb;
    PHASIC::Process_Base *proc;
    ATOOLS::Particle_Vector m_particles,  dipoleISR, dipoleFSR;
    ATOOLS::ParticleMomMap m_inparticles, m_outparticles;

    // std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR;
    ATOOLS::Poincare m_cmsboost;
    ATOOLS::Vec4D   m_Wm, m_Wp;
    Vec4D_Vector m_plab, m_reallab, m_FSRPhotons, m_ISRPhotons, m_bornMomenta, m_beams;
    ATOOLS::Flavour_Vector m_flavs;
    std::vector<double> m_mass;
    ATOOLS::Momenta_Stretcher m_stretcher;

public:

    YFS::ISR *p_isr;
    YFS::Debug *p_debug;
    ATOOLS::nlo_type::code m_nlotype;
    BEAM::Beam_Base *p_beams[2];
    YFS::Define_Dipoles *p_dipole;
    Vec4D m_photonSumISR, m_photonSumFSR;
    Vec4D_Vector m_fsrphotonsforME, m_isrphotonsforME;
    Flavour_Vector m_inflav;
    Flavour_Vector m_outflav;
    
    YFS::NLO_Base *p_nlo;
    YFS::FSR *p_fsr;
    YFS::YFS_Form_Factor *p_yfsFormFact;
    YFS::Define_Dipoles *p_dipoles;
    YFS::Ceex_Base *p_ceex;

    void Init(std::vector<double> &splimits, const Flavour_Vector &fl);

    // bool On();
    bool MakeYFS(ATOOLS::Vec4D_Vector &p);
    bool MakeYFS();
    void SetSprimeLimits(std::vector<double> &splimits);
    void SetFlavours(const Flavour_Vector &fl);
    void SetMomenta(const ATOOLS::Vec4D_Vector &plab);
    void SetBornMomenta(const ATOOLS::Vec4D_Vector &plab);
    void InitializeCEEX(const ATOOLS::Flavour_Vector &fl);
    void CalculateBeta();
    void GenerateWeight();
    void CreatMomentumMap();
    void MakeCEEX();
    Vec4D_Vector GetPhotons();
    inline double GetWeight() const      { return m_yfsweight;  }
    inline ATOOLS::Vec4D GetPhotonSum() const      { return m_photonSumISR;  }
    inline ATOOLS::Vec4D GetPhotonSumFSR() const      { return m_photonSumFSR;  }
    inline ATOOLS::Vec4D_Vector GetISRPhotons() const      { return m_ISRPhotons;  }
    inline ATOOLS::Vec4D_Vector GetFSRPhotons() const      { return m_FSRPhotons;  }
    inline ATOOLS::Flavour GetFlavour(const int i) const { return m_flavs[i];  }
    inline void SetBorn(double b) { m_born = b; }
    inline double GetBorn() { return m_born; }

    // inline void SetMode(int mode)    { m_mode = mode;}
    // inline void SetOrder(int order)  { m_order = order;}
    // inline void SetISRPhotons(Vec4D_Vector &k) {m_ISRPhotons = k;}
    inline int  GetMode()   const {return m_mode; }
    // inline int  GetBetaOrder()   const {return m_betaorder; }
    // inline int  GetOrder()   const {return m_order; }
    inline int  GetFSRMode()   const {return m_fsrmode; }
    inline int  NPhoton()   const {return m_N; }
    inline bool UsingLoopTool()   const {return m_looptool; }
    inline bool FillBlob()   const {return m_fillblob; }
    inline Vec4D_Vector GetMomenta() { return m_plab; }
    inline Vec4D GetMomenta(int i) { return m_plab[i]; }
    inline Vec4D_Vector BornMomenta() { return m_bornMomenta; }
    inline int  GetRealPhotonN()    const {return m_NReal;}
    inline void SetRealPhotonN(int n)  {m_NReal = n;}
    inline void SetQedMode(int qed)  {m_qedmode = qed;}
    inline int  GetQedMode()  {return m_qedmode;}
    inline void SetSprimeMin(double sp) {m_smin = sp;}
    inline void SetNLOType(ATOOLS::nlo_type::code nlo_type) {m_nlotype = nlo_type;}
    inline NLO_Base *NLO() {return p_nlo; }

    void SetBeam(BEAM::Beam_Spectra_Handler *beam);
    inline void SetSprime(double sp) { m_sp = sp; }
    inline void SetRunMode(const int &rmode) { m_rmode = rmode; }

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
    void YFSDebug(double W);
};// end of class YFS_Handler



}// end of namespace YFS

#endif
