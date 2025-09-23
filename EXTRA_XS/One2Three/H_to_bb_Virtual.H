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
    METOOLS::Current *m_fcur;
    std::vector<METOOLS::Vertex*> m_v1;

    std::vector<METOOLS::Current*> m_anticur;
    METOOLS::Current *m_antifcur;
    std::vector<METOOLS::Vertex*> m_antiv1;

    std::vector<size_t> m_nhel;

    PHASIC::Color_Integrator* p_ci;

    size_t NHel(const ATOOLS::Flavour& fl);
    protected:
      double alpha_qcd;
    private:
      static MODEL::Model_Base *s_model;
  public:
    H_to_bb_Virtual(const std::vector<ATOOLS::Flavour>& flavs, MODEL::Model_Base* s_model);
    ~H_to_bb_Virtual();
    void Calculate_alpha_QCD(MODEL::Model_Base* s_model);
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false) override;
    void Calculate_ME2(const std::vector<ATOOLS::Flavour>& flavs);
  };
}

