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
  class H_to_bbg_Real : public METOOLS::Spin_Amplitudes {
    size_t m_non_prop;
    size_t m_gluon;
    size_t m_propj;
    double scale;
    double alpha_qcd;

    std::vector<METOOLS::Current*> m_cur;
    METOOLS::Current *m_scur, *m_fcur;
    std::vector<METOOLS::Vertex*> m_v1, m_v2;

    std::vector<METOOLS::Current*> m_anticur;
    METOOLS::Current *m_antiscur, *m_antifcur;
    std::vector<METOOLS::Vertex*> m_antiv1, m_antiv2;

    std::vector<size_t> m_nhel;
    
    ATOOLS::Flavour m_prop;

    PHASIC::Color_Integrator* p_ci;

    size_t NHel(const ATOOLS::Flavour& fl);

    void CalculateAlphaQCD(double scale);
  public:
    H_to_bbg_Real(const std::vector<ATOOLS::Flavour>& flavs,
              const ATOOLS::Flavour& prop,
              size_t nonprop, size_t propi, size_t propj);
    ~H_to_bbg_Real();
    bool IsNLODecay();
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    double get_NLO_part();
    METOOLS::Amplitude2_Tensor AddNLOTensor(METOOLS::Amplitude2_Tensor old_tensor);
    std::string getType();
    ATOOLS::Flavour Prop() const { return m_prop; }
  };
}

