#ifndef YFS_Main_FORM_H
#define YFS_Main_FORM_H
#include "ATOOLS/Math/Vector.H"
#include "YFS/Main/YFS_Base.H"
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
    double BVR_full(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2, double omega);
    double BR_full(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2, double omega);
    double BVR_cru(double p1p2, double E1, double E2,
                   double Mas1, double Mas2, double Kmax, double MasPhot);
    double BVR_cru(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2, double E2);
    double A(double p1p2, double mass1, double mass2);
    double A(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2);
    double A4(double p1p2, double En1, double En2,
              double mass1, double mass2);
    double A4light(double p1p2, double En1, double En2,
                   double mass1, double mass2);
    double A4_eq(double E, double mass);
    double YijEta(double eta, double y1, double y2, double y3, double y4);
    double Zij(double eta, double yi, double yj);
    double Chi(double eta, double yi, double yj, double yk, double yl);
    double BVV_full(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax, int mode);
    double Full(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax, int mode);
    DivArrC BVV_full_eps(const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax, int mode);
    double BVV_WW(const ATOOLS::Vec4D_Vector born, const ATOOLS::Vec4D_Vector k, const ATOOLS::Vec4D p1, const ATOOLS::Vec4D p2, double MasPhot, double Kmax);
    double WW_t(double t, double m, double M, double k);
    double WW_s(ATOOLS::Vec4D p1, ATOOLS::Vec4D p2);
    double BVirtT(const Vec4D &p1, const Vec4D &p2 );
    double R1(const Vec4D &p1, const Vec4D &p2);
    double R2(const Vec4D &p1, const Vec4D &p2);
    // double C0(double p1, double p2, double p3);
    double C0(double p1, double p2, double p3,double m1, double m2, double m3);
    double B0(double s, double m1, double m2);
    Complex tsub(const Vec4D &p1, const Vec4D &p2, int mode, double QiQj, double theta1, double theta2);
    void WW_u();
    double  m_ww_s, m_ww_u, m_ww_t;
    double m_kmax, m_t1, m_t2, m_u1, m_u2;
    int m_useCoulomb;
    ATOOLS::Vec4D m_wp, m_wm, m_p1, m_p2;
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
