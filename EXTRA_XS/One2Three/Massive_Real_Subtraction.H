#pragma once
#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
typedef std::vector<int>    Int_Vector;

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
    ATOOLS::Vec4D_Vector born_momenta;
    double m_fullME2;

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
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    ATOOLS::Flavour Prop() const { return m_prop; }

    bool IsNLODecay();
    std::string getType();
    double get_NLO_ME2();
    void SetFullME2(double fullME2) override;
    double GetFullME2() override;
    double getColourFactor(const PHASIC::Color_Integrator* p_ci) override;
    ATOOLS::Vec4D_Vector GetMomenta() override { return born_momenta; }

    virtual void Calculate_real_subtraction(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
  
  protected:
    double ME2_Born;
    std::map<std::string, std::complex<double>> born_hel;
    ATOOLS::Flavour m_prop;
    double subtraction_term;
    double sign;

    static double V_ijk(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j,
                    ATOOLS::Vec4<double> p_k, const ATOOLS::Flavour& prop);
  };
}

