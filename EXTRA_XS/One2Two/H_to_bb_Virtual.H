#include "METOOLS/Main/Spin_Structure.H"
#include "MODEL/Main/Model_Base.H"
#include "EXTRA_XS/One2Two/NLO_Virtual.H"
#include <map>
#include <string>

namespace MODEL { 
  class Model_Base; 
}

namespace METOOLS {
  class Current;
  class Vertex;
  template<typename T> class CSpinor;
}

namespace PHASIC {
  class Color_Integrator;
}

namespace EXTRAXS {
  class H_to_bb_Virtual : public METOOLS::Spin_Amplitudes, public NLO_Virtual {
    std::vector<METOOLS::Current*> m_cur;
    std::vector<METOOLS::Current*> m_anticur;
    std::vector<size_t> m_nhel;
    static MODEL::Model_Base *s_model;
    double BornPrefactor;
    double VirtualPrefactor;
    std::array<std::complex<double>,16> M_finite;
    std::array<std::complex<double>,16> M_epsilon;
    std::array<std::complex<double>,16> M_epsilon2;
    PHASIC::Color_Integrator* p_ci;

    size_t NHel(const ATOOLS::Flavour& fl) const;
    void SetUpCurrents(const std::vector<ATOOLS::Flavour>& flavs) override;
    std::pair<std::vector<std::pair<METOOLS::CSpinor<double>*, int>>,
          std::vector<std::pair<METOOLS::CSpinor<double>*, int>>> CalculateSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti) override;
    void SetUpPrefactors(const std::vector<ATOOLS::Flavour>& flavs) override;
    void SetVirtualMatrixFinite(const ATOOLS::Vec4D_Vector& momenta) override;
    void SetVirtualMatrixE(const ATOOLS::Vec4D_Vector& momenta) override;
    void SetVirtualMatrixE2(const ATOOLS::Vec4D_Vector& momenta) override;
    std::map<std::string, std::complex<double>> CalculateV(const ATOOLS::Vec4D_Vector& momenta, bool anti=false) override;
    std::map<std::string, std::complex<double>> CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti=false) override;
  protected:
    double alpha_qcd;
  public:
    H_to_bb_Virtual(const std::vector<ATOOLS::Flavour>& flavs, MODEL::Model_Base* s_model);
    ~H_to_bb_Virtual();
    void Calculate_alpha_QCD(MODEL::Model_Base* s_model) override;
    double CalculateVirtualCorrection(const ATOOLS::Vec4D_Vector& momenta, bool anti=false) override;
    void CalculateSubtraction(const ATOOLS::Vec4D_Vector& momenta);
  };
}



