#include "METOOLS/Main/Spin_Structure.H"
#include <memory>

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
    EXTRAXS::NLO_Virtual* nlo_virtual;

    void CalculateNLO(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    size_t NHel(const ATOOLS::Flavour& fl);
    double total_virtual_correction = 0.0;
  public:
    Comix1to2(const std::vector<ATOOLS::Flavour>& flavs);
    ~Comix1to2();
    bool IsNLODecay();
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    double get_NLO_part();
  }; 

}

