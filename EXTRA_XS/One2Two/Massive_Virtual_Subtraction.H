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
  class Massive_Virtual_Subtraction : public METOOLS::Spin_Amplitudes {
    std::vector<METOOLS::Current*> m_cur;
    std::vector<METOOLS::Current*> m_anticur;
    std::vector<size_t> m_nhel;
    double BornPrefactor;
    PHASIC::Color_Integrator* p_ci;

    size_t NHel(const ATOOLS::Flavour& fl) const;
    void SetUpBornCurrents(const std::vector<ATOOLS::Flavour>& flavs);
    std::pair<std::vector<std::pair<METOOLS::CSpinor<double>*, int>>,
          std::vector<std::pair<METOOLS::CSpinor<double>*, int>>> CalculateBornSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti);
    void SetUpBornPrefactor(const std::vector<ATOOLS::Flavour>& flavs);
    std::map<std::string, std::complex<double>> CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
  protected:
    double alpha_qcd;
  public:
    Massive_Virtual_Subtraction(const std::vector<ATOOLS::Flavour>& flavs);
    ~Massive_Virtual_Subtraction();
    void CalculateAlphaQCD(double scale);
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    double CalculateFiniteSubtraction(const ATOOLS::Vec4D_Vector& momenta, double born_ME2);
    double CalculateEpsilonSubtraction(const ATOOLS::Vec4D_Vector& momenta, double born_ME2);
  };
}



