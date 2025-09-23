#include "METOOLS/Main/Spin_Structure.H"
#include "MODEL/Main/Model_Base.H"

namespace MODEL { 
  class Model_Base; 
}

namespace METOOLS {
  class Current;
  class Vertex;
}

namespace PHASIC {
  class Color_Integrator;
}

namespace EXTRAXS {
  class H_to_bb_Virtual : public METOOLS::Spin_Amplitudes {
    std::vector<METOOLS::Current*> m_cur;

    std::vector<METOOLS::Current*> m_anticur;

    std::vector<size_t> m_nhel;

    static MODEL::Model_Base *s_model;

    size_t NHel(const ATOOLS::Flavour& fl);
    void SetUpBornCurrents(const std::vector<ATOOLS::Flavour>& flavs);
    void SetUpVirtualCurrents(const std::vector<ATOOLS::Flavour>& flavs);
  protected:
    double alpha_qcd;
  public:
    H_to_bb_Virtual(const std::vector<ATOOLS::Flavour>& flavs, MODEL::Model_Base* s_model);
    ~H_to_bb_Virtual();
    void Calculate_alpha_QCD(MODEL::Model_Base* s_model);
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false) override;
    void Calculate_ME2(const std::vector<ATOOLS::Flavour>& flavs);
  };
}

