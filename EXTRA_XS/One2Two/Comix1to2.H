#include "METOOLS/Main/Spin_Structure.H"
#include <memory>
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

namespace METOOLS {
  class Current;
  class Vertex;
}

namespace PHASIC {
  class Color_Integrator;
}
namespace EXTRAXS { class H_to_bb_Virtual; }
namespace EXTRAXS { class NLO_Virtual; }

namespace EXTRAXS {

  class Comix1to2 : public METOOLS::Spin_Amplitudes {
    const std::vector<ATOOLS::Flavour>& m_flavs;
    bool isNLO;
    std::vector<METOOLS::Current*> m_cur;
    METOOLS::Current *m_fcur;
    std::vector<METOOLS::Vertex*> m_v1;

    std::vector<METOOLS::Current*> m_anticur;
    METOOLS::Current *m_antifcur;
    std::vector<METOOLS::Vertex*> m_antiv1;

    std::vector<size_t> m_nhel;

    PHASIC::Color_Integrator* p_ci;

    size_t NHel(const ATOOLS::Flavour& fl);

    // true, if a color point was already set externally
    //bool m_cpointgen;

  public:
    Comix1to2(const std::vector<ATOOLS::Flavour>& flavs);
    ~Comix1to2();
    bool IsNLODecay();
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    double get_NLO_ME2();
    std::string getType();
    METOOLS::Amplitude2_Tensor AddNLOTensor(METOOLS::Amplitude2_Tensor old_tensor);
  }; 
    //const PHASIC::Color_Integrator* GetColors() const override { return p_ci; };
    //void SetColors(const std::vector<int>& ci,
    //               const std::vector<int>& cj) override;
    //void FixColor(bool yes) override { m_cpointgen = yes; }
  //};
}

