#pragma once

#include <map>
#include <string>
#include <vector>
#include <complex>

#include "METOOLS/Main/Spin_Structure.H"
#include "MODEL/Main/Model_Base.H"

namespace METOOLS { template<typename T> class CSpinor; class Current; }
namespace ATOOLS  { class Flavour; using Vec4D_Vector = std::vector<ATOOLS::Vec4D>; }

namespace EXTRAXS {

class NLO_Virtual {
  /* this is an interface for the virtual NLO classes */
public:
  using C = std::complex<double>;
  using DDSpin = METOOLS::CSpinor<double>;
  using SpinorWithHel = std::pair<DDSpin*, int>;
  using SpinorVecWithHel = std::vector<SpinorWithHel>;
  using SpinorPairWithHel = std::pair<SpinorVecWithHel, SpinorVecWithHel>;

  virtual ~NLO_Virtual() = default;

  // setup
  virtual void Calculate_alpha_QCD(MODEL::Model_Base* s_model) = 0;
  virtual void SetUpCurrents(const std::vector<ATOOLS::Flavour>& flavs) = 0;
  virtual void SetUpPrefactors(const std::vector<ATOOLS::Flavour>& flavs) = 0;

  // spinor / matrix builders
  virtual SpinorPairWithHel CalculateSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti) = 0;
  virtual void SetVirtualMatrixFinite(const ATOOLS::Vec4D_Vector& momenta) = 0;
  virtual void SetVirtualMatrixE(const ATOOLS::Vec4D_Vector& momenta) = 0;
  virtual void SetVirtualMatrixE2(const ATOOLS::Vec4D_Vector& momenta) = 0;

  // calculation
  virtual std::map<std::string, C> CalculateV(const ATOOLS::Vec4D_Vector& momenta, bool anti=false) = 0;
  virtual std::map<std::string, C> CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti=false) = 0;
};

} // namespace EXTRAXS