#include "PDF/CJKL/CJKLph_Fortran_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Message.H"
#include <unistd.h>

using namespace PDF;
using namespace ATOOLS;

extern "C" {
void partons_(double &, double &, double *);
}

CJKLph_Fortran_Interface::CJKLph_Fortran_Interface(const ATOOLS::Flavour _bunch) {
    m_xmin = 1.e-5;
    m_xmax = 1.;
    m_q2min = 1;
    m_q2max = 2.e5;
    m_nf = 5;

    m_set = "CJKL";
    m_bunch = _bunch;
    m_d = m_u = m_s = m_c = m_b = m_g = 0.;

    for (int i = 1; i < 6; i++) {
        m_partons.insert(Flavour((kf_code)(i)));
        m_partons.insert(Flavour((kf_code)(i)).Bar());
    }
    m_partons.insert(Flavour(kf_gluon));
    m_partons.insert(Flavour(kf_jet));
    m_partons.insert(Flavour(kf_quark));
    m_partons.insert(Flavour(kf_quark).Bar());
}

PDF_Base *CJKLph_Fortran_Interface::GetCopy() {
    return new CJKLph_Fortran_Interface(m_bunch);
}

void CJKLph_Fortran_Interface::CalculateSpec(const double &_x, const double &_Q2) {
    double x = _x / m_rescale, Q2 = _Q2;

    double f[11];

    partons_(x, Q2, f);
    m_g = f[5];
    m_d = f[6];
    m_u = f[7];
    m_s = f[8];
    m_c = f[9];
    m_b = f[10];
}

double CJKLph_Fortran_Interface::GetXPDF(const ATOOLS::Flavour &infl) {
    double value = 0.;

    if (infl.Kfcode() == kf_gluon) value = m_g;
    else if (infl.Kfcode() == kf_d) value = m_d;
    else if (infl.Kfcode() == kf_u) value = m_u;
    else if (infl.Kfcode() == kf_s) value = m_s;
    else if (infl.Kfcode() == kf_c) value = m_c;
    else if (infl.Kfcode() == kf_b) value = m_b;

    value *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"),
                                            sqr(rpa->gen.Ecms()));

    return m_rescale * value;
}

double CJKLph_Fortran_Interface::GetXPDF(const kf_code &kf, bool anti) {
    double value = 0.;

    if (kf == kf_gluon) value = m_g;
    else if (kf == kf_d) value = m_d;
    else if (kf == kf_u) value = m_u;
    else if (kf == kf_s) value = m_s;
    else if (kf == kf_c) value = m_c;
    else if (kf == kf_b) value = m_b;

    value *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"),
                                            sqr(rpa->gen.Ecms()));

    return m_rescale * value;
}

DECLARE_PDF_GETTER(CJKLph_Getter);

PDF_Base *CJKLph_Getter::operator()
        (const Parameter_Type &args) const {
    // TODO: Change to another if-statement: it has to check that SAL is not specified explicitly.
    if (!args.m_bunch.IsPhoton()) return NULL;
    return new CJKLph_Fortran_Interface(args.m_bunch);
}

void CJKLph_Getter::PrintInfo
        (std::ostream &str, const size_t width) const {
    str << "CJKL photon PDF, see PRD45(1992)3986 and PRD46(1992)1973";
}

CJKLph_Getter *p_get_cjkl;

extern "C" void InitPDFLib() {
    p_get_cjkl = new CJKLph_Getter("CJKL");
}

extern "C" void ExitPDFLib() {
    delete p_get_cjkl;
}
