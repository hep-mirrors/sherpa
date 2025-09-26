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
    double BornPrefactor;
    double VirtualPrefactor;
    PHASIC::Color_Integrator* p_ci;

    size_t NHel(const ATOOLS::Flavour& fl);
    void SetUpCurrents(const std::vector<ATOOLS::Flavour>& flavs);
    void CalculateBorn(const ATOOLS::Vec4D_Vector& momenta);
    void SetUpPrefactors(const std::vector<ATOOLS::Flavour>& flavs);
    std::array<std::complex<double>,16> GetVirtualMatrixFinite(ATOOLS::Vec4D_Vector& momenta);
    std::array<std::complex<double>,16> GetVirtualMatrixE(ATOOLS::Vec4D_Vector& momenta);
    std::array<std::complex<double>,16> GetVirtualMatrixE2(ATOOLS::Vec4D_Vector& momenta);
  protected:
    double alpha_qcd;
  public:
    H_to_bb_Virtual(const std::vector<ATOOLS::Flavour>& flavs, MODEL::Model_Base* s_model);
    ~H_to_bb_Virtual();
    void Calculate_alpha_QCD(MODEL::Model_Base* s_model);
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false) override;
  };
}

