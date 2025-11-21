#include "METOOLS/Main/Spin_Structure.H"

namespace METOOLS {
  class Current;
  class Vertex;
  template<typename T> class CSpinor;
}

namespace PHASIC {
  class Color_Integrator;
}

namespace EXTRAXS {
  class H_to_bbg_Real : public METOOLS::Spin_Amplitudes {
    size_t m_non_prop;
    size_t m_gluon;
    size_t m_propj;
    double scale;
    double BornPrefactor;
    double alpha_qcd;
    std::map<std::string, std::complex<double>> born_hel;
    double ME2_Born;

    std::vector<METOOLS::Current*> m_cur;
    METOOLS::Current *m_scur, *m_fcur;
    std::vector<METOOLS::Vertex*> m_v1, m_v2;

    std::vector<METOOLS::Current*> m_anticur;
    METOOLS::Current *m_antiscur, *m_antifcur;
    std::vector<METOOLS::Vertex*> m_antiv1, m_antiv2;

    std::vector<METOOLS::Current*> m_cur_born;
    std::vector<METOOLS::Current*> m_anticur_born;

    std::vector<size_t> m_nhel;
    
    ATOOLS::Flavour m_prop;

    PHASIC::Color_Integrator* p_ci;
    PHASIC::Color_Integrator* p_ci_born;

    size_t NHel(const ATOOLS::Flavour& fl);
    double subtraction_term;

    void CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti);
    std::pair<std::vector<std::pair<METOOLS::CSpinor<double>*, int>>,
          std::vector<std::pair<METOOLS::CSpinor<double>*, int>>> CalculateSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti);
    void SetUpBornCurrents(const std::vector<ATOOLS::Flavour>& flavs);
    void SetUpBornPrefactor(const std::vector<ATOOLS::Flavour>& flavs);
    void CalculateAlphaQCD(double scale);
  public:
    H_to_bbg_Real(const std::vector<ATOOLS::Flavour>& flavs,
              const ATOOLS::Flavour& prop,
              size_t nonprop, size_t propi, size_t propj);
    ~H_to_bbg_Real();
    bool IsNLODecay();
    void Calculate_real_subtraction(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    double get_NLO_part();

    ATOOLS::Flavour Prop() const { return m_prop; }
  };
}

