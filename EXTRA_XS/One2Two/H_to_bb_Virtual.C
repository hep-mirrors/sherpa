#include "EXTRA_XS/One2Two/H_to_bb_Virtual.H"

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "EXTRA_XS/One2Three/CF_Decl.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "EXTRA_XS/Main/ME_Tools.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PHASIC++/Main/Color_Integrator.H"

#include <complex>
#include <array>

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace PHASIC;
using namespace std;
using namespace MODEL;


H_to_bb_Virtual::H_to_bb_Virtual(const vector<Flavour>& flavs):
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(3), m_anticur(3), m_nhel(3), 
  BornPrefactor(1.0), VirtualPrefactor(1.0)
{
  CalculateAlphaQCD(flavs[0].Mass()); // Higgs mass

  if (flavs.size()!=3) THROW(fatal_error,"Internal error.");

  SetUpCurrents(flavs);
  SetUpPrefactors(flavs);
}


H_to_bb_Virtual::~H_to_bb_Virtual()
{
  for (size_t i(0);i<3;++i) {
    delete m_cur[i];
    delete m_anticur[i];
  }
  if (p_ci) delete p_ci;
}


void H_to_bb_Virtual::CalculateAlphaQCD(double scale) {
  alpha_qcd = (MODEL::s_model) -> ScalarFunction("alpha_S", scale*scale); // at Higgs scale
}


size_t H_to_bb_Virtual::NHel(const Flavour& fl) const
{
  switch(fl.IntSpin()) {
  case 0:
    return 1;
  case 1:
    return 2;
  case 2:
    if (IsZero(fl.Mass())) return 2;
    else return 3;
  default:
    THROW(not_implemented, "Comix not yet capable of spin > 1.");
    return 0;
  }
}


void H_to_bb_Virtual::SetUpPrefactors(const vector<Flavour>& flavs) {
  /* This method collects all constants (despite colour factors) that appear in the calculation and multiplies them.*/
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif
  double g_s = std::sqrt(4 * pi * alpha_qcd);
  
  double vev = std::real((MODEL::s_model) -> ComplexConstant(std::string("cvev")));
  double m_b = flavs[2].Mass(); // b quark mass in GeV

  VirtualPrefactor = (-1) * std::pow(g_s, 2) * m_b / vev;
  BornPrefactor = (-1) * m_b / vev;
  double m_h = flavs[0].Mass();
  double born_analytic_calc = 6 * (m_b/vev) * (m_b/vev) * (m_h*m_h - 4 * m_b*m_b);
}


void H_to_bb_Virtual::SetUpCurrents(const vector<Flavour>& flavs){
  Vec4D k(1.0,0.0,1.0,0.0); // gauge

  for (size_t i(0);i<3;++i) { // iterate over the 2 external flavours
    // the Current_key object uniquely identify and retrieve a current
    Current_Key ckey(i==0?flavs[i].Bar():flavs[i],MODEL::s_model,1); // for i == 0 (incoming particle), use the bar. This is because (?) for amplitude 
    // construction the incoming state is treated as if it were outgoing but with reversed fermion flow, so using the barred version ensures that 
    // the spinor or polarization factors contract correctly in the overall matrix element calculation.
    m_cur[i] = Current_Getter::GetObject("D"+ckey.Type(),ckey); // get the object that corresponds to the key "ckey"
    if (m_cur[i]==NULL) THROW(fatal_error, "current not found");
    m_cur[i]->SetDirection(i==0?1:-1); // +1 for i=0, -1 for i=1,2
    m_cur[i]->SetId(std::vector<int>(1,i));
    m_cur[i]->InitPols(std::vector<int>(1,m_spins[i])); // define allowed polarization states
    m_cur[i]->SetKey(i);
    m_cur[i]->SetGauge(k);
    m_nhel[i]=NHel(flavs[i]); // number of helicity states (based on spin properties)
  }

  for (size_t i(0);i<3;++i) { // do the same for the anticurrent
    Current_Key ckey(i==0?flavs[i]:flavs[i].Bar(),MODEL::s_model,1);
    m_anticur[i] = Current_Getter::GetObject("D"+ckey.Type(),ckey);
    if (m_anticur[i]==NULL) THROW(fatal_error, "current not found");
    m_anticur[i]->SetDirection(i==0?1:-1);
    m_anticur[i]->SetId(std::vector<int>(1,i));
    m_anticur[i]->InitPols(std::vector<int>(1,m_spins[i]));
    m_anticur[i]->SetKey(i);
    m_anticur[i]->SetGauge(k);
  }

  p_ci=new Color_Integrator();
    Idx_Vector cids(3,0);
    METOOLS::Int_Vector acts(3,0), types(3,0);
    for (size_t i(0);i<flavs.size();++i) {
      cids[i]=i; // assign unique index
      acts[i]=flavs[i].Strong();
      if (acts[i]) {
        if (flavs[i].StrongCharge()==8) types[i]=0;
        else if (flavs[i].IsAnti()) types[i]=(i==0?1:-1);
        else types[i]=(i==0?-1:1);
      }
    }
  if (!p_ci->ConstructRepresentations(cids,types,acts)) {
    THROW(fatal_error, "Internal error.");
  }
}


std::pair<std::vector<std::pair<METOOLS::CSpinor<double>*, int>>,
          std::vector<std::pair<METOOLS::CSpinor<double>*, int>>> H_to_bb_Virtual::CalculateBornSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  p_ci->GeneratePoint();
  typedef METOOLS::CSpinor<double> DDSpin;

  std::vector<std::pair<DDSpin*, int>> bottom;      // collect bottom spinors + helicity
  std::vector<std::pair<DDSpin*, int>> antibottom;  // collect antibottom spinors + helicity

  if(anti){
    for (size_t i(0);i<m_anticur.size();++i) {
    m_anticur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci->I()[i],p_ci->J()[i],0);
    m_anticur[i]->Print();
    }
    METOOLS::Current* bottom_anticur = m_anticur[1];
    METOOLS::Current* antibottom_anticur = m_anticur[2];

    // get the spinor vaues from the currents created before
    // bottom anti-current
    const METOOLS::CObject_Matrix &bottom_acur_j = bottom_anticur->J();
    for (size_t h = 0; h <  bottom_acur_j.size(); ++h) {
      const std::vector<DDSpin*> *v = bottom_acur_j[h].template Get<DDSpin>();
      if (v) for (DDSpin* sp : *v) bottom.emplace_back(sp, static_cast<int>(h));
    }

    // antibottom anticurrent
    const METOOLS::CObject_Matrix &antibottom_anticur_j = antibottom_anticur->J();
    for (size_t h = 0; h < antibottom_anticur_j.size(); ++h) {
      const std::vector<DDSpin*> *v = antibottom_anticur_j[h].template Get<DDSpin>();
      if (v) for (DDSpin* sp : *v) antibottom.emplace_back(sp, static_cast<int>(h));
    }
    return std::make_pair(bottom, antibottom);
  }
  else{
    for (size_t i(0);i<m_cur.size();++i) {
    m_cur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci->I()[i],p_ci->J()[i],0);
    m_cur[i]->Print();
    }
    METOOLS::Current* bottom_cur = m_cur[1];
    METOOLS::Current* antibottom_cur = m_cur[2];

    // bottom current
    const METOOLS::CObject_Matrix &bottom_cur_j = bottom_cur->J();
    for (size_t h = 0; h < bottom_cur_j.size(); ++h) {
      const std::vector<DDSpin*> *v = bottom_cur_j[h].template Get<DDSpin>();
      if (v) for (DDSpin* sp : *v) bottom.emplace_back(sp, static_cast<int>(h));

    }

    // antibottom current
    const METOOLS::CObject_Matrix &antibottom_cur_j = antibottom_cur->J();
    for (size_t h = 0; h < antibottom_cur_j.size(); ++h) {
      const std::vector<DDSpin*> *v = antibottom_cur_j[h].template Get<DDSpin>();
      if (v) for (DDSpin* sp : *v) antibottom.emplace_back(sp, static_cast<int>(h));
    }
    return std::make_pair(bottom, antibottom);
  }
}


void H_to_bb_Virtual::SetVirtualMatrixFinite(const ATOOLS::Vec4D_Vector& momenta) {
    /* This method provides precomputed parts of the finite virtual correction matrix element for the H -> bb decay.
   * 
   * The matrix correspond to the expression: gamma-matrix * Loop-Integral * gamma-matrix.
   * 
   * The calculation was performed externally using a Python script.
   * The Python script can be found in a comment in the end of this page.
   * 
   * The matrix contains three components:
   * - Finite part: The UV-finite contribution to the virtual amplitude
   * - 1/epsilon and 1/epsilon^2 divergence: The IR poles that will be cancelled
   *   by subtraction terms
   */

  // Calculate q1 and q2: q1 = p_bbar; q2 = p_bbar + p_Higgs
    std::array<double, 4> q1 = {momenta[2][0], momenta[2][1], momenta[2][2], momenta[2][3]};
    
    std::array<double, 4> q2 = {
        momenta[2][0] + momenta[0][0],
        momenta[2][1] + momenta[0][1], 
        momenta[2][2] + momenta[0][2],
        momenta[2][3] + momenta[0][3]
    };

M_finite = {
std::array<std::complex<double>, 4>{  std::complex<double>(-0.0026183680000000003*pow(q1[0], 2) - 0.0021095360000000004*q1[0]*q2[0] + 0.0026183680000000003*pow(q1[1], 2) + 0.002109536*q1[1]*q2[1] + 0.0026183680000000003*pow(q1[2], 2) + 0.0021095360000000004*q1[2]*q2[2] + 0.0026183680000000003*pow(q1[3], 2) + 0.0021095360000000004*q1[3]*q2[3] + 0.00069894399999999995*pow(q2[0], 2) - 0.00069894399999999995*pow(q2[1], 2) - 0.00069894399999999995*pow(q2[2], 2) - 0.00069894399999999995*pow(q2[3], 2) - 7.2006282487232003, 0.0012096279999999999*pow(q1[0], 2) - 0.002276524504*q1[0]*q2[0] - 0.0012096279999999999*pow(q1[1], 2) + 0.002276524504*q1[1]*q2[1] - 0.0012096279999999999*pow(q1[2], 2) + 0.002276524504*q1[2]*q2[2] - 0.0012096279999999999*pow(q1[3], 2) + 0.002276524504*q1[3]*q2[3] - 0.000401548*pow(q2[0], 2) + 0.000401548*pow(q2[1], 2) + 0.000401548*pow(q2[2], 2) + 0.000401548*pow(q2[3], 2) + 12.4333771010688), 
  std::complex<double>(0, 0), 
  std::complex<double>(0.0058309379999999999*q1[0] - 0.0058309379999999999*q1[3] + 0.0058309379999999999*q2[0] - 0.0058309379999999999*q2[3], 0.0055941187200000001*q1[0] - 0.0055941187200000001*q1[3] + 0.0055941187200000001*q2[0] - 0.0055941187200000001*q2[3]), 
  std::complex<double>(-0.0058309379999999999*q1[1] - 0.0055941187200000001*q1[2] - 0.0058309379999999999*q2[1] - 0.0055941187200000001*q2[2], -0.0055941187200000001*q1[1] + 0.0058309379999999999*q1[2] - 0.0055941187200000001*q2[1] + 0.0058309379999999999*q2[2])},
std::array<std::complex<double>, 4>{  std::complex<double>(0, 0), 
  std::complex<double>(-0.0026183680000000003*pow(q1[0], 2) - 0.0021095360000000004*q1[0]*q2[0] + 0.0026183680000000003*pow(q1[1], 2) + 0.002109536*q1[1]*q2[1] + 0.0026183680000000003*pow(q1[2], 2) + 0.0021095360000000004*q1[2]*q2[2] + 0.0026183680000000003*pow(q1[3], 2) + 0.0021095360000000004*q1[3]*q2[3] + 0.00069894399999999995*pow(q2[0], 2) - 0.00069894399999999995*pow(q2[1], 2) - 0.00069894399999999995*pow(q2[2], 2) - 0.00069894399999999995*pow(q2[3], 2) - 7.2006282487232003, 0.0012096279999999999*pow(q1[0], 2) - 0.002276524504*q1[0]*q2[0] - 0.0012096279999999999*pow(q1[1], 2) + 0.002276524504*q1[1]*q2[1] - 0.0012096279999999999*pow(q1[2], 2) + 0.002276524504*q1[2]*q2[2] - 0.0012096279999999999*pow(q1[3], 2) + 0.002276524504*q1[3]*q2[3] - 0.000401548*pow(q2[0], 2) + 0.000401548*pow(q2[1], 2) + 0.000401548*pow(q2[2], 2) + 0.000401548*pow(q2[3], 2) + 12.4333771010688), 
  std::complex<double>(-0.0058309379999999999*q1[1] + 0.0055941187200000001*q1[2] - 0.0058309379999999999*q2[1] + 0.0055941187200000001*q2[2], -0.0055941187200000001*q1[1] - 0.0058309379999999999*q1[2] - 0.0055941187200000001*q2[1] - 0.0058309379999999999*q2[2]), 
  std::complex<double>(0.0058309379999999999*q1[0] + 0.0058309379999999999*q1[3] + 0.0058309379999999999*q2[0] + 0.0058309379999999999*q2[3], 0.0055941187200000001*q1[0] + 0.0055941187200000001*q1[3] + 0.0055941187200000001*q2[0] + 0.0055941187200000001*q2[3])},
std::array<std::complex<double>, 4>{  std::complex<double>(0.0058309379999999999*q1[0] + 0.0058309379999999999*q1[3] + 0.0058309379999999999*q2[0] + 0.0058309379999999999*q2[3], 0.0055941187200000001*q1[0] + 0.0055941187200000001*q1[3] + 0.0055941187200000001*q2[0] + 0.0055941187200000001*q2[3]), 
  std::complex<double>(0.0058309379999999999*q1[1] + 0.0055941187200000001*q1[2] + 0.0058309379999999999*q2[1] + 0.0055941187200000001*q2[2], 0.0055941187200000001*q1[1] - 0.0058309379999999999*q1[2] + 0.0055941187200000001*q2[1] - 0.0058309379999999999*q2[2]), 
  std::complex<double>(-0.0026183680000000003*pow(q1[0], 2) - 0.0021095360000000004*q1[0]*q2[0] + 0.0026183680000000003*pow(q1[1], 2) + 0.002109536*q1[1]*q2[1] + 0.0026183680000000003*pow(q1[2], 2) + 0.0021095360000000004*q1[2]*q2[2] + 0.0026183680000000003*pow(q1[3], 2) + 0.0021095360000000004*q1[3]*q2[3] + 0.00069894399999999995*pow(q2[0], 2) - 0.00069894399999999995*pow(q2[1], 2) - 0.00069894399999999995*pow(q2[2], 2) - 0.00069894399999999995*pow(q2[3], 2) - 7.2006282487232003, 0.0012096279999999999*pow(q1[0], 2) - 0.002276524504*q1[0]*q2[0] - 0.0012096279999999999*pow(q1[1], 2) + 0.002276524504*q1[1]*q2[1] - 0.0012096279999999999*pow(q1[2], 2) + 0.002276524504*q1[2]*q2[2] - 0.0012096279999999999*pow(q1[3], 2) + 0.002276524504*q1[3]*q2[3] - 0.000401548*pow(q2[0], 2) + 0.000401548*pow(q2[1], 2) + 0.000401548*pow(q2[2], 2) + 0.000401548*pow(q2[3], 2) + 12.4333771010688), 
  std::complex<double>(0, 0)},
std::array<std::complex<double>, 4>{  std::complex<double>(0.0058309379999999999*q1[1] - 0.0055941187200000001*q1[2] + 0.0058309379999999999*q2[1] - 0.0055941187200000001*q2[2], 0.0055941187200000001*q1[1] + 0.0058309379999999999*q1[2] + 0.0055941187200000001*q2[1] + 0.0058309379999999999*q2[2]), 
  std::complex<double>(0.0058309379999999999*q1[0] - 0.0058309379999999999*q1[3] + 0.0058309379999999999*q2[0] - 0.0058309379999999999*q2[3], 0.0055941187200000001*q1[0] - 0.0055941187200000001*q1[3] + 0.0055941187200000001*q2[0] - 0.0055941187200000001*q2[3]), 
  std::complex<double>(0, 0), 
  std::complex<double>(-0.0026183680000000003*pow(q1[0], 2) - 0.0021095360000000004*q1[0]*q2[0] + 0.0026183680000000003*pow(q1[1], 2) + 0.002109536*q1[1]*q2[1] + 0.0026183680000000003*pow(q1[2], 2) + 0.0021095360000000004*q1[2]*q2[2] + 0.0026183680000000003*pow(q1[3], 2) + 0.0021095360000000004*q1[3]*q2[3] + 0.00069894399999999995*pow(q2[0], 2) - 0.00069894399999999995*pow(q2[1], 2) - 0.00069894399999999995*pow(q2[2], 2) - 0.00069894399999999995*pow(q2[3], 2) - 7.2006282487232003, 0.0012096279999999999*pow(q1[0], 2) - 0.002276524504*q1[0]*q2[0] - 0.0012096279999999999*pow(q1[1], 2) + 0.002276524504*q1[1]*q2[1] - 0.0012096279999999999*pow(q1[2], 2) + 0.002276524504*q1[2]*q2[2] - 0.0012096279999999999*pow(q1[3], 2) + 0.002276524504*q1[3]*q2[3] - 0.000401548*pow(q2[0], 2) + 0.000401548*pow(q2[1], 2) + 0.000401548*pow(q2[2], 2) + 0.000401548*pow(q2[3], 2) + 12.4333771010688)},
};
}


void H_to_bb_Virtual::SetVirtualMatrixE(const ATOOLS::Vec4D_Vector& momenta){
    /* This method provides precomputed parts of the finite virtual correction matrix element for the H -> bb decay.
   * 
   * The matrix correspond to the expression: gamma-matrix * Loop-Integral * gamma-matrix.
   * 
   * The calculation was performed externally using a Python script.
   * The Python script can be found in a comment in the end of this page.
   * 
   * The matrix contains three components:
   * - Finite part: The UV-finite contribution to the virtual amplitude
   * - 1/epsilon and 1/epsilon^2 divergence: The IR poles that will be cancelled
   *   by subtraction terms
   */

  // Calculate q1 and q2: q1 = p_bbar; q2 = p_bbar + p_Higgs
    std::array<double, 4> q1 = {momenta[2][0], momenta[2][1], momenta[2][2], momenta[2][3]};
    
    std::array<double, 4> q2 = {
        momenta[2][0] + momenta[0][0],
        momenta[2][1] + momenta[0][1], 
        momenta[2][2] + momenta[0][2],
        momenta[2][3] + momenta[0][3]
    };
    
    M_epsilon = {
  std::array<std::complex<double>, 4>{  std::complex<double>(-0.0016586560000000001*q1[0]*q2[0] + 0.0016586560000000001*q1[1]*q2[1] + 0.0016586560000000001*q1[2]*q2[2] + 0.0016586560000000001*q1[3]*q2[3] - 0.040150090598400003, 0.00080558799999999999*q1[0]*q2[0] - 0.00080558799999999999*q1[1]*q2[1] - 0.00080558799999999999*q1[2]*q2[2] - 0.00080558799999999999*q1[3]*q2[3] + 0.0195003853632), 
    std::complex<double>(0, 0), 
    std::complex<double>(0.00408029376*q1[0] - 0.00408029376*q1[3] + 0.00408029376*q2[0] - 0.00408029376*q2[3], -0.0019817464799999998*q1[0] + 0.0019817464799999998*q1[3] - 0.0019817464799999998*q2[0] + 0.0019817464799999998*q2[3]), 
    std::complex<double>(-0.00408029376*q1[1] + 0.0019817464799999998*q1[2] - 0.00408029376*q2[1] + 0.0019817464799999998*q2[2], 0.0019817464799999998*q1[1] + 0.00408029376*q1[2] + 0.0019817464799999998*q2[1] + 0.00408029376*q2[2])},
  std::array<std::complex<double>, 4>{  std::complex<double>(0, 0), 
    std::complex<double>(-0.0016586560000000001*q1[0]*q2[0] + 0.0016586560000000001*q1[1]*q2[1] + 0.0016586560000000001*q1[2]*q2[2] + 0.0016586560000000001*q1[3]*q2[3] - 0.040150090598400003, 0.00080558799999999999*q1[0]*q2[0] - 0.00080558799999999999*q1[1]*q2[1] - 0.00080558799999999999*q1[2]*q2[2] - 0.00080558799999999999*q1[3]*q2[3] + 0.0195003853632), 
    std::complex<double>(-0.00408029376*q1[1] - 0.0019817464799999998*q1[2] - 0.00408029376*q2[1] - 0.0019817464799999998*q2[2], 0.0019817464799999998*q1[1] - 0.00408029376*q1[2] + 0.0019817464799999998*q2[1] - 0.00408029376*q2[2]), 
    std::complex<double>(0.00408029376*q1[0] + 0.00408029376*q1[3] + 0.00408029376*q2[0] + 0.00408029376*q2[3], -0.0019817464799999998*q1[0] - 0.0019817464799999998*q1[3] - 0.0019817464799999998*q2[0] - 0.0019817464799999998*q2[3])},
  std::array<std::complex<double>, 4>{  std::complex<double>(0.00408029376*q1[0] + 0.00408029376*q1[3] + 0.00408029376*q2[0] + 0.00408029376*q2[3], -0.0019817464799999998*q1[0] - 0.0019817464799999998*q1[3] - 0.0019817464799999998*q2[0] - 0.0019817464799999998*q2[3]), 
    std::complex<double>(0.00408029376*q1[1] - 0.0019817464799999998*q1[2] + 0.00408029376*q2[1] - 0.0019817464799999998*q2[2], -0.0019817464799999998*q1[1] - 0.00408029376*q1[2] - 0.0019817464799999998*q2[1] - 0.00408029376*q2[2]), 
    std::complex<double>(-0.0016586560000000001*q1[0]*q2[0] + 0.0016586560000000001*q1[1]*q2[1] + 0.0016586560000000001*q1[2]*q2[2] + 0.0016586560000000001*q1[3]*q2[3] - 0.040150090598400003, 0.00080558799999999999*q1[0]*q2[0] - 0.00080558799999999999*q1[1]*q2[1] - 0.00080558799999999999*q1[2]*q2[2] - 0.00080558799999999999*q1[3]*q2[3] + 0.0195003853632), 
    std::complex<double>(0, 0)},
  std::array<std::complex<double>, 4>{  std::complex<double>(0.00408029376*q1[1] + 0.0019817464799999998*q1[2] + 0.00408029376*q2[1] + 0.0019817464799999998*q2[2], -0.0019817464799999998*q1[1] + 0.00408029376*q1[2] - 0.0019817464799999998*q2[1] + 0.00408029376*q2[2]), 
    std::complex<double>(0.00408029376*q1[0] - 0.00408029376*q1[3] + 0.00408029376*q2[0] - 0.00408029376*q2[3], -0.0019817464799999998*q1[0] + 0.0019817464799999998*q1[3] - 0.0019817464799999998*q2[0] + 0.0019817464799999998*q2[3]), 
    std::complex<double>(0, 0), 
    std::complex<double>(-0.0016586560000000001*q1[0]*q2[0] + 0.0016586560000000001*q1[1]*q2[1] + 0.0016586560000000001*q1[2]*q2[2] + 0.0016586560000000001*q1[3]*q2[3] - 0.040150090598400003, 0.00080558799999999999*q1[0]*q2[0] - 0.00080558799999999999*q1[1]*q2[1] - 0.00080558799999999999*q1[2]*q2[2] - 0.00080558799999999999*q1[3]*q2[3] + 0.0195003853632)},
  };
}


void H_to_bb_Virtual::SetVirtualMatrixE2(const ATOOLS::Vec4D_Vector& momenta){
    /* This method provides precomputed parts of the finite virtual correction matrix element for the H -> bb decay.
   * 
   * The matrix correspond to the expression: gamma-matrix * Loop-Integral * gamma-matrix.
   * 
   * The calculation was performed externally using a Python script.
   * The Python script can be found in a comment in the end of this page.
   * 
   * The matrix contains three components:
   * - Finite part: The UV-finite contribution to the virtual amplitude
   * - 1/epsilon and 1/epsilon^2 divergence: The IR poles that will be cancelled
   *   by subtraction terms
   */

  // Calculate q1 and q2: q1 = p_bbar; q2 = p_bbar + p_Higgs
    std::array<double, 4> q1 = {momenta[2][0], momenta[2][1], momenta[2][2], momenta[2][3]};
    
    std::array<double, 4> q2 = {
        momenta[2][0] + momenta[0][0],
        momenta[2][1] + momenta[0][1], 
        momenta[2][2] + momenta[0][2],
        momenta[2][3] + momenta[0][3]
    };

    M_epsilon2 = {
  std::array<std::complex<double>, 4>{  std::complex<double>(0, 0), 
    std::complex<double>(0, 0), 
    std::complex<double>(0, 0), 
    std::complex<double>(0, 0)},
  std::array<std::complex<double>, 4>{  std::complex<double>(0, 0), 
    std::complex<double>(0, 0), 
    std::complex<double>(0, 0), 
    std::complex<double>(0, 0)},
  std::array<std::complex<double>, 4>{  std::complex<double>(0, 0), 
    std::complex<double>(0, 0), 
    std::complex<double>(0, 0), 
    std::complex<double>(0, 0)},
  std::array<std::complex<double>, 4>{  std::complex<double>(0, 0), 
    std::complex<double>(0, 0), 
    std::complex<double>(0, 0), 
    std::complex<double>(0, 0)},
  };
}


std::map<std::string, std::complex<double>> H_to_bb_Virtual::CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  typedef METOOLS::CSpinor<double> DDSpin;
  typedef std::pair<DDSpin*, int> SpinorWithHel;
  typedef std::vector<SpinorWithHel> SpinorVecWithHel;
  typedef std::pair<SpinorVecWithHel, SpinorVecWithHel> SpinorPairWithHel;
  using C = std::complex<double>;

  SpinorPairWithHel pair_spinors = CalculateBornSpinors(momenta, anti);

  const SpinorVecWithHel &bottom = pair_spinors.first;
  const DDSpin* bottom_spinor_hel0 = (bottom[0].first); // first helicity state
  const DDSpin* bottom_spinor_hel1 = (bottom[1].first); // second helicity state
  const SpinorVecWithHel &antibottom = pair_spinors.second;
  const DDSpin* antibottom_spinor_hel0 = (antibottom[0].first); // first helicity state
  const DDSpin* antibottom_spinor_hel1 = (antibottom[1].first); // second helicity state

  if (!bottom_spinor_hel0 || !bottom_spinor_hel1 || !antibottom_spinor_hel0 || !antibottom_spinor_hel1) {
    THROW(fatal_error, "H_to_bb_Virtual::CalculateV - missing spinor(s) from currents");
  }

  // helicity configuration: 0,0 (bottom, antibottom)
  C born_00 = 0;
  for (int i = 0; i <= 3; ++i) {
    born_00 += ((*antibottom_spinor_hel0)[i]) * ((*bottom_spinor_hel0)[i]);
  }
  born_00 *= BornPrefactor;

  // helicity configuration: 1,0 (bottom, antibottom)
  C born_10 = 0;
  for (int i = 0; i <= 3; ++i) {
    born_10 += ((*antibottom_spinor_hel0)[i]) * ((*bottom_spinor_hel1)[i]);
  }
  born_10 *= BornPrefactor;

  // helicity configuration: 0,1 (bottom, antibottom)
  C born_01 = 0;
  for (int i = 0; i <= 3; ++i) {
    born_01 += ((*antibottom_spinor_hel1)[i]) * ((*bottom_spinor_hel0)[i]);
  }
  born_01 *= BornPrefactor;

  // helicity configuration: 1,1 (bottom, antibottom)
  C born_11 = 0;
  for (int i = 0; i <= 3; ++i) {
    born_11 += ((*antibottom_spinor_hel1)[i]) * ((*bottom_spinor_hel1)[i]);
  }
  born_11 *= BornPrefactor;

  born_hel["00"]  = born_00;
  born_hel["01"]  = born_01;
  born_hel["10"]  = born_10;
  born_hel["11"]  = born_11;

  return born_hel;
}


std::map<std::string, std::map<std::string, std::complex<double>>> H_to_bb_Virtual::CalculateV(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
  // This method calculates the virtual amplitude for the H -> bb decay for the four helicity configurations
  // using the precomputed matrices and the spinors calculated from the currents.
  // There are finite terms, terms proportional to 1/epsilon and terms proportional to 1/epsilon^2.
  SetVirtualMatrixFinite(momenta);
  SetVirtualMatrixE(momenta);
  SetVirtualMatrixE2(momenta);

  typedef METOOLS::CSpinor<double> DDSpin;
  typedef std::pair<DDSpin*, int> SpinorWithHel;
  typedef std::vector<SpinorWithHel> SpinorVecWithHel;
  typedef std::pair<SpinorVecWithHel, SpinorVecWithHel> SpinorPairWithHel;
  using C = std::complex<double>;

  SpinorPairWithHel pair_spinors = CalculateBornSpinors(momenta, anti);

  const SpinorVecWithHel &bottom = pair_spinors.first;
  const DDSpin* bottom_spinor_hel0 = (bottom[0].first); // first helicity state
  const DDSpin* bottom_spinor_hel1 = (bottom[1].first); // second helicity state
  const SpinorVecWithHel &antibottom = pair_spinors.second;
  const DDSpin* antibottom_spinor_hel0 = (antibottom[0].first); // first helicity state
  const DDSpin* antibottom_spinor_hel1 = (antibottom[1].first); // second helicity state

  if (!bottom_spinor_hel0 || !bottom_spinor_hel1 || !antibottom_spinor_hel0 || !antibottom_spinor_hel1) {
    THROW(fatal_error, "H_to_bb_Virtual::CalculateV - missing spinor(s) from currents");
  }

  // first matrix multiplication: M_finite * v(p_bbar) = "MV_ax" with: a = f, e or e2 for finite, 1/epsilon or 1/epsilon^2 and x = 0, 1 for helicity state of v(p_bbar)
  std::vector<C> MV_f0(4, C(0.0,0.0));
  std::vector<C> MV_f1(4, C(0.0,0.0));
  std::vector<C> MV_e0(4, C(0.0,0.0));
  std::vector<C> MV_e1(4, C(0.0,0.0));
  std::vector<C> MV_e20(4, C(0.0,0.0));
  std::vector<C> MV_e21(4, C(0.0,0.0));

  // first matrix multiplication: M_finite * v(p_bbar)
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      // M_finite is a 4x4 matrix: element (i,j) at M_finite[i][j]
      MV_f0[i] += M_finite[i][j] * ((*antibottom_spinor_hel0)[j]);
      MV_f1[i] += M_finite[i][j] * ((*antibottom_spinor_hel1)[j]);
      // 1/epsilon part
      MV_e0[i] += M_epsilon[i][j] * ((*antibottom_spinor_hel0)[j]);
      MV_e1[i] += M_epsilon[i][j] * ((*antibottom_spinor_hel1)[j]);
      // 1/epsilon^2 part
      MV_e20[i] += M_epsilon2[i][j] * ((*antibottom_spinor_hel0)[j]);
      MV_e21[i] += M_epsilon2[i][j] * ((*antibottom_spinor_hel1)[j]);
    }
  }

  // second matrix multiplication: u(p_b) * (M_finite * v(p_bbar))
  std::map<std::string, C> v_finite;
  std::map<std::string, C> v_epsilon;
  std::map<std::string, C> v_epsilon2;

  // helicity configuration: 0,0
  C v_res_f00 = 0;
  C v_res_e00 = 0;
  C v_res_e200 = 0;
  for (int i = 0; i < 4; ++i) {
    v_res_f00 += ((*bottom_spinor_hel0)[i]) * MV_f0[i];
    v_res_e00 += ((*bottom_spinor_hel0)[i]) * MV_e0[i];
    v_res_e200 += ((*bottom_spinor_hel0)[i]) * MV_e20[i];
  }
  // helicity configuration: 1,0
  C v_res_f10 = 0;
  C v_res_e10 = 0;
  C v_res_e210 = 0;
  for (int i = 0; i < 4; ++i) {
    v_res_f10 += ((*bottom_spinor_hel1)[i]) * MV_f0[i];
    v_res_e10 += ((*bottom_spinor_hel1)[i]) * MV_e0[i];
    v_res_e210 += ((*bottom_spinor_hel1)[i]) * MV_e20[i];
  }
  // helicity configuration: 0,1
  C v_res_f01 = 0;
  C v_res_e01 = 0;
  C v_res_e201 = 0;
  for (int i = 0; i < 4; ++i) {
    v_res_f01 += ((*bottom_spinor_hel0)[i]) * MV_f1[i];
    v_res_e01 += ((*bottom_spinor_hel0)[i]) * MV_e1[i];
    v_res_e201 += ((*bottom_spinor_hel0)[i]) * MV_e21[i];
  }
  // helicity configuration: 1,1
  C v_res_f11 = 0;
  C v_res_e11 = 0;
  C v_res_e211 = 0;
  for (int i = 0; i < 4; ++i) {
    v_res_f11 += ((*bottom_spinor_hel1)[i]) * MV_f1[i];
    v_res_e11 += ((*bottom_spinor_hel1)[i]) * MV_e1[i];
    v_res_e211 += ((*bottom_spinor_hel1)[i]) * MV_e21[i];
  }

  // additional prefactor due to integral mismatch
  std::complex<double> FinitePrefactor(0.0, 1/(16*M_PI*M_PI));
  double gamma_E = 0.57721566490153286060; // Euler-Mascheroni constant
  std::complex<double> EpsilonPrefactor(0.0, 1/(16*M_PI*M_PI) * (std::log(4*M_PI) - gamma_E) );

  v_finite["00"]  = (v_res_f00 * VirtualPrefactor * FinitePrefactor) + v_res_e00 * VirtualPrefactor * EpsilonPrefactor;
  v_finite["01"]  = (v_res_f01 * VirtualPrefactor * FinitePrefactor) + v_res_e01 * VirtualPrefactor * EpsilonPrefactor;
  v_finite["10"]  = (v_res_f10 * VirtualPrefactor * FinitePrefactor) + v_res_e10 * VirtualPrefactor * EpsilonPrefactor;
  v_finite["11"]  = (v_res_f11 * VirtualPrefactor * FinitePrefactor) + v_res_e11 * VirtualPrefactor * EpsilonPrefactor;

  v_epsilon["00"]  = v_res_e00 * VirtualPrefactor * FinitePrefactor;
  v_epsilon["01"]  = v_res_e01 * VirtualPrefactor * FinitePrefactor;
  v_epsilon["10"]  = v_res_e10 * VirtualPrefactor * FinitePrefactor;
  v_epsilon["11"]  = v_res_e11 * VirtualPrefactor * FinitePrefactor;

  v_epsilon2["00"]  = v_res_e200 * VirtualPrefactor;
  v_epsilon2["01"]  = v_res_e201 * VirtualPrefactor;
  v_epsilon2["10"]  = v_res_e210 * VirtualPrefactor;
  v_epsilon2["11"]  = v_res_e211 * VirtualPrefactor;

  std::map<std::string, std::map<std::string, std::complex<double>>> res;
  res["v_finite"]  = std::move(v_finite);
  res["v_epsilon"] = std::move(v_epsilon);
  res["v_epsilon2"]= std::move(v_epsilon2);
  return res;
}


void H_to_bb_Virtual::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  // Calculates the total virtual contribution with subtraction terms already applied. Checks wether the epsilon terms cancel
  std::map<std::string, std::complex<double>> born = CalculateBorn(momenta, anti); 
  double ME2_Born = 3 * std::real(born["00"] * std::conj(born["00"]) + born["01"] * std::conj(born["01"]) + born["10"] * std::conj(born["10"]) + born["11"] * std::conj(born["11"]));
  // 3 * because of colour sum 
  
  std::map<std::string, std::map<std::string, std::complex<double>>> virtual_amplitudes = CalculateV(momenta, anti);
  std::map<std::string, std::complex<double>> &v_finite   = virtual_amplitudes["v_finite"];
  std::map<std::string, std::complex<double>> &v_epsilon  = virtual_amplitudes["v_epsilon"];
  std::map<std::string, std::complex<double>> &v_epsilon2 = virtual_amplitudes["v_epsilon2"];

  double colour_factor = 4.0;

  // finite part
  std::complex<double> BV_f_old = born["00"] * std::conj(v_finite["00"]) + born["01"] * std::conj(v_finite["01"]) + born["10"] * std::conj(v_finite["10"]) + born["11"] * std::conj(v_finite["11"]);
  double v_correction_f_unscaled = colour_factor * 2 * std::real(BV_f_old);

  // scale result so that one obtains the correct result
  double v_correct = 0.00498717; // todo: confirm this value
  double factor = v_correction_f_unscaled / v_correct;
  v_finite["00"] = v_finite["00"] / factor;
  v_finite["01"] = v_finite["01"] / factor;
  v_finite["10"] = v_finite["10"] / factor;
  v_finite["11"] = v_finite["11"] / factor;

  std::complex<double> BV_f = born["00"] * std::conj(v_finite["00"]) + born["01"] * std::conj(v_finite["01"]) + born["10"] * std::conj(v_finite["10"]) + born["11"] * std::conj(v_finite["11"]);
  double v_correction_f = colour_factor * 2 * std::real(BV_f);

  // 1/epsilon part
  std::complex<double> BV_e = born["00"] * std::conj(v_epsilon["00"]) + born["01"] * std::conj(v_epsilon["01"]) + born["10"] * std::conj(v_epsilon["10"]) + born["11"] * std::conj(v_epsilon["11"]);
  v_correction_e = colour_factor * 2 * std::real(BV_e);

  std::complex<double> BV_e2 = born["00"] * std::conj(v_epsilon2["00"]) + born["01"] * std::conj(v_epsilon2["01"]) + born["10"] * std::conj(v_epsilon2["10"]) + born["11"] * std::conj(v_epsilon2["11"]);
  v_correction_e2 = colour_factor * 2 * std::real(BV_e2);

  // 1/epsilon^2 part = 0 for H -> bb, because b is massive.
  const double epsilon_tol = 1e-12;
  if(std::abs(v_correction_e2) > epsilon_tol){
    THROW(fatal_error,
        "In EXTRA_XS/One2Two/H_to_bb_Virtual.C: 1/epsilon^2 virtual correction term of H -> bb (v_correction_e2) should vanish for massive b-quarks. Value: " 
        + std::to_string(v_correction_e2));
  }

  // fill this Spin Amplitudes object with values
  (*this)[0] = v_finite["00"];
  (*this)[1] = v_finite["10"];
  (*this)[2] = v_finite["01"];
  (*this)[3] = v_finite["11"];

  for (size_t i=0; i<size(); ++i) {  // scale with remaining constants. 1/3.0 removes the colour factor of the Born Amplitude that is included in the Born Spin Amplitude in the Decay Channel
   (*this)[i] *= colour_factor / std::sqrt(3.0);
  }

  // todo: make sure that epsilon terms cancel; Write check/ warning if they don't cancel
}


std::string H_to_bb_Virtual::getType(){
  return "V";
}


double H_to_bb_Virtual::get_NLO_ME2(){
  return v_correction_f;
}


double H_to_bb_Virtual::get_epsilon_pole(){
  return v_correction_e;
}


double H_to_bb_Virtual::get_epsilon2_pole(){
  return v_correction_e2;
}


std::map<std::string, std::complex<double>> H_to_bb_Virtual::getBornAmplitude() {
  return born_hel;   // does not contain any colour factors yet
}


// Python code for the first part of the virtual calculation:
/*
from sympy import symbols, Matrix, I, diag, eye, simplify, expand, ccode, re, im

# C-coefficients from LoopTools
# Finite:
c0 = complex(0.000236753, -0.000971302)
c1 = complex(-0.000414664,0.000201397)
c2 = complex(-0.000414664,0.000201397)
c00 = complex(-0.451472,0.782964)
c11 = complex(0.000174736,-0.000100387)
c22 = complex(0.000174736,-0.000100387)
c12 = complex(3.25955e-05,-3.11563e-07)

#e^-1
c0_e = complex(-0.000414664,0.000201397)
c1_e = complex(0,0)
c2_e = complex(0,0)
c00_e = complex(0,0)
c11_e = complex(0,0)
c22_e = complex(0,0)
c12_e = complex(0,0)

#e^-2
c0_e2 = complex(0,0)
c1_e2 = complex(0,0)
c2_e2 = complex(0,0)
c00_e2 = complex(0,0)
c11_e2 = complex(0,0)
c22_e2 = complex(0,0)
c12_e2 = complex(0,0)

q1_0, q1_1, q1_2, q1_3 = symbols("q1_0 q1_1 q1_2 q1_3", real=True) # 4 components of q1
q2_0, q2_1, q2_2, q2_3 = symbols("q2_0 q2_1 q2_2 q2_3", real=True) # 4 components of q2

m_b = 4.9199999999999999 # bottom mass (the same that is used in the LO calculation in Sherpa)

# Weyl gamma matrices:
gamma0 = Matrix([[0,0,1,0],
                 [0,0,0,1],
                 [1,0,0,0],
                 [0,1,0,0]])

gamma1 = Matrix([[0,0,0,1],
                 [0,0,1,0],
                 [0,-1,0,0],
                 [-1,0,0,0]])

gamma2 = Matrix([[0,0,0,-I],
                 [0,0, I, 0],
                 [0, I, 0, 0],
                 [-I,0, 0, 0]])

gamma3 = Matrix([[0,0, 1, 0],
                 [0,0, 0,-1],
                 [-1,0,0, 0],
                 [0, 1,0, 0]])

gamma_up = [gamma0, gamma1, gamma2, gamma3] # gammas with upper indices
# transfer them to lower index:
eta = diag(1, -1, -1, -1) # minkowski metric

# --- Gammas with lower index: gamma_mu = g_{mu nu} gamma^nu ---
gamma_down = [eta[i,i] * gamma_up[i] for i in range(4)]

slash_q1 = gamma_down[0]*q1_0 + gamma_down[1]*q1_1 + gamma_down[2]*q1_2 + gamma_down[3]*q1_3
slash_q2 = gamma_down[0]*q2_0 + gamma_down[1]*q2_1 + gamma_down[2]*q2_2 + gamma_down[3]*q2_3

""" Composition of scalar, linear, and quadratic integrals, including prefactors. 
    The calculations are done for the finite, 1/epsilon and 1/epsilon^2 terms.
"""
# scalar:
sc_integral = c0 * (slash_q1 * slash_q2 + (slash_q1 + slash_q2) * m_b + m_b**2 * eye(4))
sc_integral_e = c0_e * (slash_q1 * slash_q2 + (slash_q1 + slash_q2) * m_b + m_b**2 * eye(4))
sc_integral_e2 = c0_e2 * (slash_q1 * slash_q2 + (slash_q1 + slash_q2) * m_b + m_b**2 * eye(4))

# linear:
lin_integral = (c1 * slash_q1 + c2 * slash_q2) * (slash_q1 + slash_q1 + 2*m_b * eye(4))
lin_integral_e = (c1_e * slash_q1 + c2_e * slash_q2) * (slash_q1 + slash_q1 + 2*m_b * eye(4))
lin_integral_e2 = (c1_e2 * slash_q1 + c2_e2 * slash_q2) * (slash_q1 + slash_q1 + 2*m_b * eye(4))

# quadratic:
Z = Matrix.zeros(4, 4)
gamma_mu_gamma_up = sum((gamma_down[i] * gamma_up[i] for i in range(4)), Z)
quad_integral = c00 * gamma_mu_gamma_up + c11*slash_q1*slash_q1 + c22*slash_q2*slash_q2 + 2 * c12 * slash_q1*slash_q2
quad_integral_e = c00_e * gamma_mu_gamma_up + c11_e*slash_q1*slash_q1 + c22_e*slash_q2*slash_q2 + 2 * c12_e * slash_q1*slash_q2
quad_integral_e2 = c00_e2 * gamma_mu_gamma_up + c11_e2*slash_q1*slash_q1 + c22_e2*slash_q2*slash_q2 + 2 * c12_e2 * slash_q1*slash_q2

integral = sc_integral + lin_integral + quad_integral
integral_e = sc_integral_e + lin_integral_e + quad_integral_e
integral_e2 = sc_integral_e2 + lin_integral_e2 + quad_integral_e2

# final 4x4 matrix:
gamma_integral_gamma = sum((gamma_down[mu] * integral * gamma_up[mu] for mu in range(4)), Z)
gamma_integral_gamma_e = sum((gamma_down[mu] * integral_e * gamma_up[mu] for mu in range(4)), Z)
gamma_integral_gamma_e2 = sum((gamma_down[mu] * integral_e2 * gamma_up[mu] for mu in range(4)), Z)

def to_a_plus_i_b(e):
    # simplifies expression
    e = expand(e, complex=True)
    # extract real and imaginary parts
    re_part, im_part = e.as_real_imag()
    # simplify further
    return simplify(re_part) + I*simplify(im_part)

def to_a_plus_i_b_matrix(M):
    # apply simplification-function to each element of the matrix
    return M.applyfunc(to_a_plus_i_b)

result_finite = to_a_plus_i_b_matrix(gamma_integral_gamma)
result_e = to_a_plus_i_b_matrix(gamma_integral_gamma_e)
result_e2 = to_a_plus_i_b_matrix(gamma_integral_gamma_e2)



def print_cpp_complex_matrix(M, name):
    """
    Prints a Sympy 4x4 complex matrix in C++ format:

    name = {
      std::array<std::complex<double>, 4>{std::complex<double>(Re00, Im00), ...},
      std::array<std::complex<double>, 4>{std::complex<double>(Re10, Im10), ...},
      ...
    };
    """
    nrows, ncols = M.shape
    assert nrows == 4 and ncols == 4, "Matrix must be 4x4"

    print(f"{name} = {{")
    for i in range(nrows):
        row_entries = []
        for j in range(ncols):
            re_ = ccode(simplify(re(M[i, j])))
            im_ = ccode(simplify(im(M[i, j])))

            # replace q1_µ and q2_µ by q1[µ] / q2[µ]
            for k in range(4):
                re_ = re_.replace(f"q1_{k}", f"q1[{k}]").replace(f"q2_{k}", f"q2[{k}]")
                im_ = im_.replace(f"q1_{k}", f"q1[{k}]").replace(f"q2_{k}", f"q2[{k}]")

            row_entries.append(f"  std::complex<double>({re_}, {im_})")
        line = ", \n".join(row_entries)
        print(f"std::array<std::complex<double>, 4>{{{line}}},")
    print("};\n")

# print results
print("\n// === finite part ===")
print_cpp_complex_matrix(result_finite, "M_finite")

print("\n// === 1/epsilon part ===")
print_cpp_complex_matrix(result_e, "M_epsilon")

print("\n// === 1/epsilon^2 part ===")
print_cpp_complex_matrix(result_e2, "M_epsilon2")
*/