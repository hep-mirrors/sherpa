#ifndef YFS_Main_FORM_H
#define YFS_Main_FORM_H
#include "ATOOLS/Math/Vector.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/Dipole.H"
#include "YFS/Main/Leg.H"
#include "YFS/Main/Soft_Region.H"
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
    //! The region-explicit forms (Soft_Region.H): leg energies are taken in the
    //! region's frame. The energy-based overloads above serve YFS/Tools checks.
    double BVR_full(const SoftRegion &R, const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2,
                    double Mas1, double Mas2, double MasPhot, int mode);
    double BVR_cru(const SoftRegion &R, const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2,
                   double Mas1, double Mas2);
    //! real B-tilde over R plus the dipole's virtual (the FormFactorSum term)
    double BVR_full(YFS::Dipole &d, const SoftRegion &R);
    //! the exponent region of a II/FF dipole at cutoff omega: a final pair in
    //! its own rest frame (FF_FORMFACTOR_REST_FRAME), everything else in P
    SoftRegion ExponentRegion(const YFS::Dipole &d, double omega) const;
    double BVR_decay(double M, double ml, double ks);

    DivArrD BVR_full_eps(YFS::Dipole &d,  double Kmax, int mode);
    double BVR_full(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2, double omega);
    double BVR_cru(double p1p2, const Leg &l1, const Leg &l2, double Kmax);
    double BVR_cru(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2, double E2);
    double A(double p1p2, double mass1, double mass2);
    double A(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2);
    double A1(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2);
    double A2(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2);
    double A4(double p1p2, Leg l1, Leg l2);
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
    double BVirtT(Vec4D p1, Vec4D p2,  double kmax=0);
    //! BVirtT in dim. reg., with BVV_full_eps's map of ln(lambda^2): the virtual B of an
    //! initial-final (t-channel) dipole, the one the exponent's IFForFac carries
    DivArrD BVirtT_eps(YFS::Dipole &d);
    // YFS exponent for an initial-final dipole: Btilda plus the t-channel
    // virtual (KKMC's TForFac).
    double IFForFac(YFS::Dipole &d, double omega);
    double IFForFac(YFS::Dipole &d, const SoftRegion &R);
    double  BVirtGeneral(YFS::Dipole &d, double kmax=0);
    DivArrD BVirtGeneralEps(YFS::Dipole &d, double kmax=0);
    // Virtual YFS B from the dim-reg master integrals (B0, C0) rather than the
    // analytic A1/A2 forms; carries the 1/eps_IR pole (FULL_FORM==3).
    DivArrD BVirtMaster(YFS::Dipole &d, double kmax=0);
    double  m_ww_s, m_ww_u, m_ww_t;
    double m_t1, m_t2, m_u1, m_u2;
    int m_useCoulomb;
    ATOOLS::Vec4D m_wp, m_wm, m_p1, m_p2;
    //! Not owned: the loop provider (set by NLO_Base), read for its IR scale
    //! and epsilon scheme. Null when unset.
    PHASIC::Virtual_ME2_Base *p_virt = nullptr;
    double m_m1, m_m2;
    //! The couplings for a pair of charge product ch (FSR::Initialize, once
    //! per dipole and event), set relative to the unit-charge values.
    inline void SetCharge(double ch) {
        if (!m_unitcharge) { m_alpha_unit = m_alpha; m_alpi_unit = m_alpi; m_unitcharge = true; }
        m_alpha = m_alpha_unit * ch;
        m_alpi = m_alpi_unit * ch;
    }
    bool   m_unitcharge = false;
    double m_alpha_unit = 0., m_alpi_unit = 0.;
};

}


#endif
