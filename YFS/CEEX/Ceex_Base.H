#ifndef YFS_CEEX_Ceex_handler
#define YFS_CEEX_Ceex_handler

#include "ATOOLS/Phys/Flavour.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "EXTAMP/Process.H"
#include "YFS/Main/YFS_Base.H"


using namespace MODEL;
using namespace ATOOLS;
using namespace std;
using namespace METOOLS;

namespace YFS {

  typedef std::map<int,double> HelicityMap;
  typedef complex<double> dcmplx;

  class Amplitude{
  public:
    Amplitude();
    ~Amplitude(){}
    Complex m_A[2][2][2][2];
    Complex m_U[2][2];
    Complex m_V[2][2];
  };

  class Ceex_Base : public YFS::YFS_Base {  
  public:

    Ceex_Base(const Flavour_Vector &flavs);

    Flavour_Vector m_flavs;
    PHASIC::Virtual_ME2_Base*     p_loop_me;
    PHASIC::Tree_ME2_Base*        p_lo_me;
    Amplitude  m_bornAmp, m_Amp1, m_Amp1U, m_Amp1V;
    Vec4D_Vector   m_momenta, m_isrphotons,m_bornmomenta;
    Vec4D m_eta, m_zeta, m_b, m_pp1, m_pp2, m_kk;
    int m_onlyz, m_onlyg;
    bool m_fixedwidth,m_checkxs;
    std::vector<int> m_PhoHel;
    PHASIC::Process_Info m_pi, m_pireal;
    XYZFunc    *p_xyz;
    Complex    m_T,m_U, m_sW, m_cW, m_sin2tw,m_gcpl, m_zcpl, m_propG, m_propZ, m_bornsum, m_cfac;
    Complex m_cL, m_cR, m_prop, m_I_L, m_I_R, m_F_L, m_F_R, m_I, m_Soft;
    Complex m_Tamp[2][2][2][2], m_Uamp[2][2][2][2],m_ampborn[2][2][2][2], m_Sc, m_Tc, m_Uc;
    Complex m_Sini[2][2], m_TC[2], m_UC[2], m_ampISR1,  m_ae, m_af, m_ve, m_vf, m_norm;
    Complex m_BoxGGtu, m_BoxGZtu, m_BoxGGut, m_BoxGZut, m_Sprod, m_beta00, m_beta10, m_beta01, m_beta20;
    // Complex m_Amp1U[2][2], m_Amp1V[2][2];
    Complex m_vertexI, m_b1, m_b2;
    double     m_Q1Q2I, m_Q1Q2F, m_QIQF, m_MZ, m_gZ, m_isrformfactor, m_sQ;
    double     m_mass_I, m_mass_F, m_qe, m_qf, m_crude;
    double     m_result, m_e;
    void RegisterDefaults();
    void Init(const Vec4D_Vector &p);
    virtual void Calculate();
    void InfraredSubtractedME_0_0();
    void InfraredSubtractedME_0_1();
    void InfraredSubtractedME_1_0(Vec4D &k, int hel);
    void InfraredSubtractedME_2_0();
    Complex BetaDouble_2_0(Vec4D &k1, Vec4D &k2, int h1, int h2);
    Complex BetaSingle_2_0(Vec4D &k1, Vec4D &k2, int h1, int h2);
    Complex BetaRest_2_0(Vec4D &k1, Vec4D &k2, int h1, int h2);
    void CalculateSfactors();
    void MakePhotonHel();
    void LoadME();
    Complex BoxGG();
    Complex BoxGZ();
    Complex BoxSubtract();
    Complex Soft(Vec4D k, Vec4D p1, Vec4D p2, int hel1);
    Complex Sfactor(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int hel);
    Complex UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h2, int i, int j, bool NegMass=0);
    Complex VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h2, int i, int j);
    Complex bsigma(Vec4D k, Vec4D p, int hel);
    virtual inline bool ZOnly() {return m_onlyz;}

    Complex S(const Vec4D &p1, const Vec4D &p2, int h1, int h2);
    Complex S(const Vec4D &p1, const Vec4D &p2, double m1,double m2, int h1, int h2);
    Complex T(const Vec4D &p1, const Vec4D &p2, int h1, int h2);
    Complex Tp(const Vec4D &p1, const Vec4D &p2, int h1, int h2);
    Complex U(const Vec4D &p1, const Vec4D &p2, int h1, int h2);
    Complex Up(const Vec4D &p1, const Vec4D &p2, int h1, int h2);
    Complex T_mass(const Vec4D &p1, const Vec4D &p2, double m1, double m2, int h1, int h2);
    Complex Tp_mass(const Vec4D &p1, const Vec4D &p2, double m1, double m2, int h1, int h2);
    Complex U_mass(const Vec4D &p1, const Vec4D &p2,  double m1, double m2, int h1, int h2);
    Complex Up_mass(const Vec4D &p1, const Vec4D &p2, double m1, double m2, int h1, int h2);
    Complex BornAmplitude(const Vec4D_Vector &k);
    Complex BornAmplitude_mass(const Vec4D_Vector &k, double m1, double m2, double m3, double m4);
    Complex BornAmplitude(const Vec4D_Vector &k, int h0, int h1, int h2, int h3);
    Complex BornAmplitude_mass(const Vec4D_Vector &k,double m1, double m2, double m3, double m4, int h0, int h1, int h2, int h3);
    Complex BornAmplitude(Vec4D p1, Vec4D p2, Vec4D p3, Vec4D p4, int h0, int h1, int h2, int h3);
    void BornAmplitude(const Vec4D_Vector &k, Amplitude &T);
    int MapHel(int &h);
    void UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpU);
    void VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpV);

    void SumAmplitude(Complex &sum, const Amplitude &Amp, const Complex fac=1.);
    void AddU(Complex &sum, const Amplitude &Amp, const Amplitude &U, const  Complex fac=1.);
    void AddV(Complex &sum, const Amplitude &Amp, const Amplitude &V, const  Complex fac=1.);
    void SumAmplitude(Complex &sum, const Amplitude &Amp1, const Amplitude &Amp2, const  Complex fac=1.);

    Complex Splus(const Vec4D &p1, const Vec4D &p2);
    Complex Sminus(const Vec4D &p1, const Vec4D &p2);
    Complex CouplingZ(double i, int mode);
    Complex CouplingG();
    void MakeProp();
    void Reset();
    
    inline double GetResult() {return m_result;}
    inline void SetISRPhotons(Vec4D_Vector &k) {m_isrphotons=k;}
    inline void SetBornMomenta (Vec4D_Vector &p) {m_bornmomenta =p;}
    inline void SetBorn (const double &born) {m_born=born;}
    inline void SetISRFormFactor (double &form) {m_isrformfactor =form;}
    inline void SetProcessInfo(PHASIC::Process_Info  const &pi) {m_pi = pi;}
    inline bool CheckXs() {return m_checkxs;}
    double Xi(const Vec4D p, const Vec4D q);
    Complex iProd(const int i, const Vec4D &p, const Vec4D &q);
  };
}

#endif
