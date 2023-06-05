#ifndef YFS_NLO_REAL_FF_H
#define YFS_NLO_REAL_FF_H
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Math/Histogram.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/Tools/Dipole.H"
#include "YFS/Main/Define_Dipoles.H"

using namespace ATOOLS;
// #define YFS_DEBUG_REAL
namespace YFS {class XS_ee_ffbar;}
namespace PHASIC {class Tree_ME2_Base;}

// namespace EXTRAXS {class ME2_Base;}
namespace YFS{
#ifdef YFS_DEBUG_REAL
  typedef std::map<std::string,ATOOLS::Histogram *> histograms;
  static histograms m_histograms;
#endif



  class Real_ff : public YFS_Base {
  public:
    Real_ff(int order);//, const PHASIC::External_ME_Args args);
    ~Real_ff();
    Vec4D m_p1, m_p2, m_q1, m_q2, m_beam1, m_beam2, m_photonSum;
    Vec4D_Vector m_photons;
    // std::map<ATOOLS::Vec4D, double> m_BETA20, m_BETA21;
    // EXTRAXS::ME2_Base *p_ME;
    PHASIC::Tree_ME2_Base *p_tree;
    int m_order, m_fsrmode, m_usefulleik, m_formfactor, m_usenllv, m_use_fac;
    double m_alp, m_alpi, m_L, m_gamma, m_gammap, m_Emin, m_u;
    double m_mass, m_mass2, m_p1p2, m_cth11, m_cth12, m_cth21, m_cth22;
    double m_born1, m_born2,m_borncrude, m_D, m_D10, m_D11, m_D12, m_beta1;
    double m_born11, m_born12, m_born21,m_born22;
    double m_beta00, m_beta01, m_beta02, m_beta03, m_beta01f, m_beta01i, m_beta02f, m_beta02i;
    double  m_beta20, m_beta20Int, m_beta21, m_beta22, m_beta23, m_yBeta10;
    double m_beta30, m_spp;
    double  m_gammaI, m_gammaF, m_delI2, m_delF2;
    double m_D20, m_D21, m_D22, m_D30;
    double m_real, m_virtual, m_realISR, m_realFSR;
    double m_vlim1, m_vlim2, m_vlim3, m_g1, m_g2, m_gg1,m_gg2, m_ggg1, m_ggg2;
    bool m_bornInit = false, m_setvirt=false, m_use_model_alpha, m_use_nnlo;
    std::map<std::pair<int,int>, double> m_BETA20, m_BETA21;
    std::vector<double> m_sisr, m_sfsr;
    std::vector<double> m_a;
    std::vector<double> m_b;
    std::vector<double> m_dvir;
    std::vector<double> m_y, m_yisr, m_zisr, m_yfsr, m_zfsr;
    std::vector<double> m_z;
    std::vector<double> m_beta10, m_beta10i, m_beta10f;
    std::vector<double> m_beta11, m_beta12, m_beta11i, m_beta11f, m_beta12i, m_beta12f, m_beta21i, m_beta21f;
    std::vector<double> m_Sfac, m_Hfac, m_d10vec;
    void SetIncoming(YFS::Dipole_Vector::iterator dipole,Vec4D_Vector &beams,Vec4D_Vector &k,std::vector<double> y, std::vector<double> z);
    void SetIncoming(YFS::Dipole *d);
    double Eikonal(Vec4D k);
    double W(double a, double b);
    double Eik(Vec4D p1, Vec4D p2, Vec4D k);
    double chi(double u, double a, double b);
    double chi3(double u, double a1, double b1, double a2, double b2);
    void Beta1(Vec4D k, double a, double b);
    void Beta2();
    double Beta10(Vec4D k1);
    double Beta10(int i);
    double Beta11(Vec4D k1, int i);
    double Beta12(Vec4D k1, int i);
    double Beta11(Vec4D k1);
    double Beta12(Vec4D k1);
    // double Beta20(Vec4D k1, Vec4D k2);
    double Beta20(Vec4D k1, Vec4D k2, int i, int k);
    double IntIF(Vec4D_Vector &k1, Vec4D_Vector &k2, std::vector<double> yi, std::vector<double> zi, std::vector<double> yf, std::vector<double> zf);
    double Beta20Int(Vec4D k1, Vec4D k2, int i, int j,  int mode);
    // double Beta21(Vec4D k1, Vec4D k2);
    double Beta21(Vec4D k1, Vec4D k2, int i, int j);
    void Beta3();
    void Calculate();
    void   D1(Vec4D k,double a, double b, double wm, int order); //Eq 14 in https://arxiv.org/pdf/hep-ph/0006359.pdf
    void   D2(int i, int j); //Eq 14 in https://arxiv.org/pdf/hep-ph/0006359.pdf
    void D3(int i, int j,  int k);
    void D3Int(int i, int j,  int k, int mode);
    double Delta(int order, double z);
    double Alpha(Vec4D kj);
    double Beta(Vec4D ki);
    void  SetVirtual(std::vector<double> vir);
    void  Calc(int n);
    double Calc(const ATOOLS::Vec4D_Vector &p);
    void GetBeta02();
    inline double GetBeta1() {return m_beta1;}
    inline double GetReal() {return m_real;}
    inline void SetBorn( double born ) {m_born = m_beta00 = born;}
    inline void SetMode( int mode ) {m_fsrmode = mode;}
    inline void SetGammaI( double g ) {m_gammaI = g;}
    inline void SetGammaF( double g ) {m_gammaF = g;}
    inline void SetSprime( double s ) {m_sp = s;}
    inline double Virtual() { return m_realISR*m_realFSR; }
    inline Vec4D GetPhotons() {return m_photonSum;}
    inline void SetSpp(double spp) { m_spp = spp;}
    inline void SetBeta01(double v) {m_beta01 = v;}
    void Sort(Vec4D_Vector &p);
    double operator()(const ATOOLS::Vec4D_Vector& momenta);
    void ResetReal();
    double CalculateBeta(const Vec4D& p);
    void CalculateVirt();
    double AddVirtual();
    double AddVirtual(int a);
    double EikonalInterferance(Vec4D k);
    double wm0(double x, double y);
    double wmd(double x, double y);
    double wmd(double del, double x, double y);
    // inline PHASIC::Tree_ME2_Base\* SetTree( PHASIC::Tree_ME2_Base* T) {p_tree = T;}
    // bool SetColours(const ATOOLS::Vec4D_Vector& momenta);
  };
}






#endif
