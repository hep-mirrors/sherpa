#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

namespace METOOLS {
  class Current;
  class Vertex;
}

namespace PHASIC {
  class Color_Integrator;
}

namespace EXTRAXS {
  class Comix1to3 : public METOOLS::Spin_Amplitudes {
    std::vector<METOOLS::Current*> m_cur;
    METOOLS::Current *m_scur, *m_fcur;
    std::vector<METOOLS::Vertex*> m_v1, m_v2;

    std::vector<METOOLS::Current*> m_anticur;
    METOOLS::Current *m_antiscur, *m_antifcur;
    std::vector<METOOLS::Vertex*> m_antiv1, m_antiv2;

    std::vector<size_t> m_nhel;
    
    ATOOLS::Flavour m_prop;

    PHASIC::Color_Integrator* p_ci;
    ATOOLS::Vec4D_Vector m_momenta;

    size_t NHel(const ATOOLS::Flavour& fl);

    // true, if a color point was already set externally
    bool m_cpointgen;
  public:
    Comix1to3(const std::vector<ATOOLS::Flavour>& flavs,
              const ATOOLS::Flavour& prop,
              size_t nonprop, size_t propi, size_t propj);
    ~Comix1to3();
    bool IsNLODecay();
    double get_NLO_ME2();
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    std::string getType();
    ATOOLS::Flavour Prop() const { return m_prop; }
    void SetColors(const std::vector<int>& ci,
                   const std::vector<int>& cj) override;
    void FixColor(bool yes) override { m_cpointgen = yes; }
    ATOOLS::Vec4D_Vector GetMomenta() override { return m_momenta; }
  };
}

