#ifndef YFS_Main_FORM_H
#define YFS_Main_FORM_H
#include "ATOOLS/Math/Vector.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/Dipole.H"
#include "YFS/NLO/Virtual.H"
#include "METOOLS/Loops/Divergence_Array.H"

using namespace METOOLS;

namespace YFS {
class YFS_Form_Factor : public YFS_Base
{
public:
    YFS_Form_Factor();
    ~YFS_Form_Factor();
    double BVR_full(double p1p2, double E1, double E2,
                    double Mas1, double Mas2, double Kmax, double MasPhot, int mode);
    double BVR_full(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2,  double Kmax, double MasPhot, int mode);
    double BVR_full(YFS::Dipole &d, double omega);
    DivArrD BVR_full_eps(YFS::Dipole &d,  double Kmax, int mode);
    double BVR_full(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2, double omega);
    double BVR_cru(double p1p2, double E1, double E2,
                   double Mas1, double Mas2, double Kmax);
    double BVR_cru(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2, double E2);
    double A(double p1p2, double mass1, double mass2);
    double A(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2);
    double A1(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2);
    double A2(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2);
    double A4(double p1p2, double En1, double En2,
              double mass1, double mass2);
    double A4_eq(double E, double mass);
    double YijEta(double eta, double y1, double y2, double y3, double y4);
    double Zij(double eta, double yi, double yj);
    double Chi(double eta, double yi, double yj, double yk, double yl);
    double BVV_full(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax, int mode);
    double BVV_full(YFS::Dipole &d, double MasPhot, double Kmax, int mode);
    // shared body of the two BVV_full overloads; masses passed in explicitly
    double BVV_full_impl(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2,
                         double Mas1, double Mas2, double MasPhot, int mode);
    DivArrD BVV_full_eps(YFS::Dipole &d, double Kmax, int mode);
    double BVV_WW(const ATOOLS::Vec4D_Vector born, const ATOOLS::Vec4D_Vector k, const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax);
    double WW_t(double t, double m, double M, double k);
    double WW_s(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2);
    double BVirtT(YFS::Dipole &d, double kmax=0);
    DivArrD BVirtTEps(YFS::Dipole &d, double kmax=0);
    double BVirtT(Vec4D p1, Vec4D p2,  double kmax=0);
    // YFS exponent for an INITIAL-FINAL (t-channel-like) dipole: Btilda plus the
    // T-CHANNEL virtual. This is KKMC's TForFac (SRCee/KKceex.cxx), which builds
    // its interference form factor Yint out of exactly exp(Btilda + TBvirt).
    double IFForFac(YFS::Dipole &d, double omega);
    double R1(YFS::Dipole &d);
    double  BVirtGeneral(YFS::Dipole &d, double kmax=0);
    DivArrD BVirtGeneralEps(YFS::Dipole &d, double kmax=0);
    // Virtual YFS B factor built from the dim-reg master integrals
    // (METOOLS/Loops: Master_Bubble B0 + Master_Triangle C0) instead of the
    // analytic A1/A2 dilog forms. Returns the DivArr (carrying the 1/eps_IR
    // pole) so it cancels against the real BVR_full_eps in the FULL_FORM==3
    // form-factor mode.
    DivArrD BVirtMaster(YFS::Dipole &d, double kmax=0);
    void WW_u();
    double  m_ww_s, m_ww_u, m_ww_t;
    double m_kmax, m_t1, m_t2, m_u1, m_u2;
    int m_useCoulomb;
    ATOOLS::Vec4D m_wp, m_wm, m_p1, m_p2;
    PHASIC::Virtual_ME2_Base *p_virt;
    double m_xx1, m_xx2, m_m1, m_m2;
    bool m_use_model_alpha;
    inline void SetCharge(double ch) {
        m_alpha = m_alpha * ch;
        m_alpi = m_alpi * ch;
    }
    inline void SetCoulomb(int c) {m_useCoulomb = c;}
};

}


#endif
