#pragma once
#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

namespace METOOLS {
  class Current;
  class Vertex;
  template<typename T> class CSpinor;
}

namespace PHASIC {
  class Color_Integrator;
}

namespace EXTRAXS {
  class Massive_Real_Subtraction : public METOOLS::Spin_Amplitudes {
    size_t m_non_prop;
    size_t m_gluon;
    size_t m_propj;
    double scale;
    double BornPrefactor;
    double alpha_qcd;

    std::vector<METOOLS::Current*> m_cur_born;
    std::vector<METOOLS::Current*> m_anticur_born;

    std::vector<size_t> m_nhel;

    PHASIC::Color_Integrator* p_ci_born;

    size_t NHel(const ATOOLS::Flavour& fl);

    void CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti);
    std::pair<std::vector<std::pair<METOOLS::CSpinor<double>*, int>>,
          std::vector<std::pair<METOOLS::CSpinor<double>*, int>>> CalculateBornSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti);
    void SetUpBornCurrents(const std::vector<ATOOLS::Flavour>& flavs);
    void SetUpBornPrefactor(const std::vector<ATOOLS::Flavour>& flavs);
    void CalculateAlphaQCD(double scale);
  public:
    Massive_Real_Subtraction(const std::vector<ATOOLS::Flavour>& flavs,
              const ATOOLS::Flavour& prop,
              size_t nonprop, size_t propi, size_t propj);
    ~Massive_Real_Subtraction();
    bool IsNLODecay();
    virtual void Calculate_real_subtraction(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    double get_NLO_ME2();
    std::string getType();
    ATOOLS::Flavour Prop() const { return m_prop; }
  
  protected:
    double ME2_Born;
    std::map<std::string, std::complex<double>> born_hel;
    ATOOLS::Flavour m_prop;
    double subtraction_term;

    static double V_ijk(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j,
                    ATOOLS::Vec4<double> p_k, const ATOOLS::Flavour& prop);
  };
}

