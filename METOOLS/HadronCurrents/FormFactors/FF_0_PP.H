#ifndef METOOLS_Current_Library_FF_0_PP_H
#define METOOLS_Current_Library_FF_0_PP_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  enum FF_0_PP_mode {
    pipi_plus       = 1,
    KK_plus         = 2,
    etapi_plus      = 3,
    etaprimepi_plus = 4,
    Kpi_plus        = 11,
    Keta_plus       = 12,
    Ketaprime_plus  = 13,
    pipi_0          = 101,
    KK_0            = 102,
    etapi_0         = 103,
    etaprimepi_0    = 104,
    unknown         = 999
  };
  
  class FF_0_PP_Base : public FormFactor_Base {
  protected:
    FF_0_PP_mode        m_mode;
    double              m_norm;
    // True if the ACTUAL final-state kaon is K_S or K_L rather than
    // the bare (unobservable) K^0/K0bar flavour eigenstate. K_S/K_L
    // are strangeness-eigenstate MIXTURES, K_{S,L}=(K^0 -+ K0bar)/
    // sqrt(2); since the weak current here only ever produces one
    // specific flavour content (K0bar, never K^0, for Cabibbo-
    // suppressed Kpi/KK - cf. the analogous discussion for the 3-meson
    // K_S pi K_S etc. currents), projecting onto K_S or K_L each pick
    // up an ADDITIONAL 1/sqrt(2) amplitude factor on top of whatever
    // isospin Clebsch-Gordan is already in m_norm for the "K0bar+X"
    // vs "K-+X'" charge assignment - set by FixMode(), applied in
    // FixParameters() of the derived classes.
    bool                 m_isKSKL;
    Summed_Propagator * p_props;
    // Family-fallback propagator (request #2): built from the plain
    // "family base" model (KS=100 for the 100s, RChiPT=200 for the
    // 200s) whenever the SPECIFICALLY requested model has no
    // dedicated construction for this channel. Used as a second tier
    // in operator(), before the final constant-only fallback. NULL
    // for ff_model::KS/RChiPT themselves (nothing to fall back to
    // that isn't p_props itself), and for channels where even the
    // family base has no implementation.
    Summed_Propagator * p_props_fallback;

    void    FixMode();
    virtual Complex FF_RChiPT(const double & Q2);
    virtual Complex FF_RChL2012(const double & Q2) { return Complex(0.,0.); }
  public:
    FF_0_PP_Base(const FF_Parameters & params) :
      FormFactor_Base(params),
      m_mode(FF_0_PP_mode::unknown), m_isKSKL(false),
      m_norm(0.), p_props(NULL), p_props_fallback(NULL) {
      FixMode();
    }
    ~FF_0_PP_Base() {
      if (p_props)          delete p_props;
      if (p_props_fallback)  delete p_props_fallback;
    }

    Complex operator()(const ATOOLS::Vec4D_Vector& moms);
    inline const double & GetNorm() const { return m_norm; }
  };
}

#endif
