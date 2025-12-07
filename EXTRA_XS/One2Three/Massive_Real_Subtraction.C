#include "EXTRA_XS/One2Three/Massive_Real_Subtraction.H"

#include "METOOLS/Explicit/Current.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include <assert.h>
#include "ATOOLS/Math/Vec4.H"
#include "EXTRA_XS/One2Three/CF_Decl.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace PHASIC;
using namespace std;

Massive_Real_Subtraction::Massive_Real_Subtraction(const vector<Flavour>& flavs, const Flavour& prop,
                     size_t non_prop, size_t gluon, size_t propj) :
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur_born(3), m_anticur_born(3), m_nhel(4), m_prop(prop), m_non_prop(non_prop), m_gluon(gluon), m_propj(propj)
{
  std::vector<Flavour> born_flavs;
  for (size_t i = 0; i < flavs.size(); ++i) {
    if(flavs[i].IDName() == "G"){
      continue;
    }
    else{
      born_flavs.push_back(flavs[i]);
    }
} 
  scale = flavs[0].Mass(); // scale = mass of incoming particle
  CalculateAlphaQCD(scale);
  SetUpBornCurrents(born_flavs);
  SetUpBornPrefactor(born_flavs);
}


Massive_Real_Subtraction::~Massive_Real_Subtraction()
{
  for (size_t i(0);i<3;++i) {
    delete m_cur_born[i];
    delete m_anticur_born[i];
  }
  if (p_ci_born) delete p_ci_born;
}


void Massive_Real_Subtraction::SetUpBornCurrents(const vector<Flavour>& flavs){
  Vec4D k(1.0,0.0,1.0,0.0); // gauge

  for (size_t i(0);i<3;++i) { // iterate over the 2 external flavours
    // the Current_key object uniquely identify and retrieve a current
    Current_Key ckey(i==0?flavs[i].Bar():flavs[i],MODEL::s_model,1); // for i == 0 (incoming particle), use the bar. This is because (?) for amplitude 
    // construction the incoming state is treated as if it were outgoing but with reversed fermion flow, so using the barred version ensures that 
    // the spinor or polarization factors contract correctly in the overall matrix element calculation.
    m_cur_born[i] = Current_Getter::GetObject("D"+ckey.Type(),ckey); // get the object that corresponds to the key "ckey"
    if (m_cur_born[i]==NULL) THROW(fatal_error, "current not found");
    m_cur_born[i]->SetDirection(i==0?1:-1); // +1 for i=0, -1 for i=1,2
    m_cur_born[i]->SetId(std::vector<int>(1,i));
    m_cur_born[i]->InitPols(std::vector<int>(1,m_spins[i])); // define allowed polarization states
    m_cur_born[i]->SetKey(i);
    m_cur_born[i]->SetGauge(k);
    m_nhel[i]=NHel(flavs[i]); // number of helicity states (based on spin properties)
  }

  for (size_t i(0);i<3;++i) { // do the same for the anticurrent
    Current_Key ckey(i==0?flavs[i]:flavs[i].Bar(),MODEL::s_model,1);
    m_anticur_born[i] = Current_Getter::GetObject("D"+ckey.Type(),ckey);
    if (m_anticur_born[i]==NULL) THROW(fatal_error, "current not found");
    m_anticur_born[i]->SetDirection(i==0?1:-1);
    m_anticur_born[i]->SetId(std::vector<int>(1,i));
    m_anticur_born[i]->InitPols(std::vector<int>(1,m_spins[i]));
    m_anticur_born[i]->SetKey(i);
    m_anticur_born[i]->SetGauge(k);
  }

  p_ci_born=new Color_Integrator();
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
  if (!p_ci_born->ConstructRepresentations(cids,types,acts)) {
    THROW(fatal_error, "Internal error.");
  }
}


void Massive_Real_Subtraction::SetUpBornPrefactor(const vector<Flavour>& flavs) {
  /* This method collects all constants (despite colour factors) that appear in the calculation and multiplies them.*/
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif
  double g_s = std::sqrt(4 * pi * alpha_qcd);
  double G_F = 1.16637886e-5; // Fermi constant in GeV^-2; value from PDG
  double vev = 1 / std::sqrt(G_F * std::sqrt(2)); // vacuum expectation value in GeV; doublecheck that value
  double m_b = flavs[2].Mass(); // b quark mass in GeV

  BornPrefactor = (-1) * m_b / vev;
}


void Massive_Real_Subtraction::CalculateAlphaQCD(double scale) {
  alpha_qcd = (MODEL::s_model) -> ScalarFunction("alpha_S", scale*scale); // at Higgs scale
}


bool Massive_Real_Subtraction::IsNLODecay(){
  return true;
}


void Massive_Real_Subtraction::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
  DEBUG_FUNC(momenta.size());
  CalculateBorn(momenta, anti);
  Calculate_real_subtraction(momenta, anti);
}


// in the following there are some helper functions defined for the real subtraction
static double v_pq(ATOOLS::Vec4<double> p, ATOOLS::Vec4<double> q){
  double pq = p * q;
  return std::sqrt(1-(p*p * q*q)/(pq * pq));
}


static double z_i_tilde(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  return p_i * p_k / (p_i * p_k + p_j * p_k);
}


static double z_j_tilde(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  return p_j * p_k / (p_j * p_k + p_i * p_k);
}


static double y_ijk(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  return p_i * p_j / (p_i * p_j + p_i * p_k + p_j * p_k);
}


static double mu_n(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k, double m_n){
  ATOOLS::Vec4<double> Q = p_i + p_j+ p_k;
  return m_n / (std::sqrt(Q*Q));
}


static double nu_ijk(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  double m_j = std::sqrt(p_j * p_j);
  double m_k = std::sqrt(p_k * p_k);
  double mu_j = mu_n(p_i, p_j, p_k, m_j);
  double mu_k = mu_n(p_i, p_j, p_k, m_k);
  double bracket = 2 * mu_k * mu_k + (1 - mu_j*mu_j - mu_k*mu_k)*(1 - y_ijk(p_i, p_j, p_k));
  double numerator = std::sqrt(bracket * bracket - 4 * mu_k * mu_k);
  double denominator = (1 - mu_j*mu_j - mu_k*mu_k)*(1 - y_ijk(p_i, p_j, p_k));
  return numerator/denominator;
}


static double lambda(double x, double y, double z){
  return x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
}


static double nu_ijk_tilde(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  double m_i = std::sqrt(p_i * p_i);
  double m_j = std::sqrt(p_j * p_j);
  double m_k = std::sqrt(p_k * p_k);
  // either i or j is the gluon and massless. Therefore, the mass m_ij equals the bottom mass. Check, which mass is the bottom mass:
  double m_ij;
  if((p_i * p_i) < 0.00001){
    m_ij = m_j; // mass so small that p_i is the gluon
  }
  else{
    m_ij = m_i;
  }
  double mu_ij = mu_n(p_i, p_j, p_k, m_ij);
  double mu_k = mu_n(p_i, p_j, p_k, m_k);
  double lambda_val = lambda(1, mu_ij*mu_ij, mu_k*mu_k);
  return std::sqrt(lambda_val)/ (1 - mu_ij*mu_ij - mu_k*mu_k);
}


static double V_ijk(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j,
                    ATOOLS::Vec4<double> p_k, const ATOOLS::Flavour& prop){
  double C_F = 4.0 / 3.0;
  ATOOLS::Vec4<double> Q = p_i + p_j + p_k;
  double m2_Q = Q*Q;
  double m2_quark = sqr(prop.Mass()); // Q mass squared
  double alpha_qcd = MODEL::s_model -> ScalarFunction("alpha_S", m2_Q); // at Higgs scale
  
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif

  double prefactor = 8 * pi * alpha_qcd * C_F;
  double first_summand = 2 / (1 - z_j_tilde(p_i, p_j, p_k) * (1 - y_ijk(p_i, p_j, p_k) ));
  double second_summand = nu_ijk_tilde(p_i, p_j, p_k) / nu_ijk(p_i, p_j, p_k) * (1 + z_j_tilde(p_i, p_j, p_k) + m2_quark / (p_i * p_j));

  return prefactor * (first_summand - second_summand);
}


void Massive_Real_Subtraction::Calculate_real_subtraction(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  // implementation based on the formulas in the Catani Dittmaier Seymour Trocsanyi paper from 2002

  // first: some variables needed later (nomination also based on paper)
  ATOOLS::Vec4<double> p_g = momenta[1];
  ATOOLS::Vec4<double> p_b = momenta[2];
  ATOOLS::Vec4<double> p_bb = momenta[3];

  double V_gb_bb = V_ijk(p_g, p_b, p_bb, m_prop);
  double V_gbb_b = V_ijk(p_g, p_bb, p_b, m_prop);

  double m2_ij = p_b * p_b; // because m_i = 0 (gluon) and m_b = m_bb

  double D_gb_bb = V_gb_bb/ ((p_g + p_b)*(p_g + p_b) - m2_ij) * ME2_Born;
  double D_gbb_b = V_gbb_b/ ((p_g + p_bb)*(p_g + p_bb) - m2_ij) * ME2_Born;

  subtraction_term = D_gb_bb + D_gbb_b;
}


double Massive_Real_Subtraction::get_NLO_part(){
  return (-1) * subtraction_term;
}


size_t Massive_Real_Subtraction::NHel(const Flavour& fl)
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


static ATOOLS::Vec4<double> p_ij_tilde(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  // mapping of momenta from real to subtraction kinematics
  ATOOLS::Vec4<double> Q = p_i + p_j + p_k;
  double Q2 = Q*Q;

  double m_i = std::sqrt(p_i * p_i);
  double m_j = std::sqrt(p_j * p_j);
  double m_ij;
  if((p_i * p_i) < 0.00001){
    m_ij = m_j; // mass so small that p_i is the gluon
  }
  else{
    m_ij = m_i;
  }
  double m_k = std::sqrt(p_k*p_k);

  double factor1 = std::sqrt(lambda(Q2, m_ij*m_ij, m_k*m_k)) / std::sqrt(lambda(Q2, (p_i + p_j) * (p_i + p_j), m_k*m_k));
  double factor2 = (Q2 + m_k*m_k - m_ij*m_ij) / (2 * Q2);
  return factor1 * (p_k - (Q*p_k/Q2) * Q) + factor2 * Q;
}


// calculate the Born Matrix Element:
void Massive_Real_Subtraction::CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  typedef METOOLS::CSpinor<double> DDSpin;
  typedef std::pair<DDSpin*, int> SpinorWithHel;
  typedef std::vector<SpinorWithHel> SpinorVecWithHel;
  typedef std::pair<SpinorVecWithHel, SpinorVecWithHel> SpinorPairWithHel;
  using C = std::complex<double>;

  ATOOLS::Vec4D_Vector born_momenta(3);
  born_momenta[0] = momenta[0];                       // incoming momentum stays the same
  // map the final momenta to the emitter and spectator momenta (formulas taken from Catani Seymour paper)
  ATOOLS::Vec4<double> Q = momenta[m_gluon] + momenta[m_propj] + momenta[m_non_prop];
  ATOOLS::Vec4<double> var_p_ij_tilde = p_ij_tilde(momenta[m_gluon], momenta[m_propj], momenta[m_non_prop]); // emitter
  ATOOLS::Vec4<double> p_k_tilde = Q - var_p_ij_tilde; // spectator
  born_momenta[1] = var_p_ij_tilde;  
  born_momenta[2] = p_k_tilde;                

  SpinorPairWithHel pair_spinors = CalculateBornSpinors(born_momenta, anti);

  const SpinorVecWithHel &bottom = pair_spinors.first;
  const DDSpin* bottom_spinor_hel0 = (bottom[0].first); // first helicity state
  const DDSpin* bottom_spinor_hel1 = (bottom[1].first); // second helicity state
  const SpinorVecWithHel &antibottom = pair_spinors.second;
  const DDSpin* antibottom_spinor_hel0 = (antibottom[0].first); // first helicity state
  const DDSpin* antibottom_spinor_hel1 = (antibottom[1].first); // second helicity state

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

  ME2_Born = 3 * std::real(born_hel["00"] * std::conj(born_hel["00"]) + born_hel["01"] * std::conj(born_hel["01"]) + born_hel["10"] * std::conj(born_hel["10"]) + born_hel["11"] * std::conj(born_hel["11"]));
  // 3 * because of colour sum
}


std::pair<std::vector<std::pair<METOOLS::CSpinor<double>*, int>>,
          std::vector<std::pair<METOOLS::CSpinor<double>*, int>>> Massive_Real_Subtraction::CalculateBornSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  p_ci_born->GeneratePoint();
  typedef METOOLS::CSpinor<double> DDSpin;

  std::vector<std::pair<DDSpin*, int>> bottom;      // collect bottom spinors + helicity
  std::vector<std::pair<DDSpin*, int>> antibottom;  // collect antibottom spinors + helicity

  if(anti){
    for (size_t i(0);i<m_anticur_born.size();++i) {
    m_anticur_born[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci_born->I()[i],p_ci_born->J()[i],0);
    m_anticur_born[i]->Print();
    }
    METOOLS::Current* bottom_anticur = m_anticur_born[1];
    METOOLS::Current* antibottom_anticur = m_anticur_born[2];

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
    for (size_t i(0);i<m_cur_born.size();++i) {
    m_cur_born[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci_born->I()[i],p_ci_born->J()[i],0);
    m_cur_born[i]->Print();
    }
    METOOLS::Current* bottom_cur = m_cur_born[1];
    METOOLS::Current* antibottom_cur = m_cur_born[2];

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


METOOLS::Amplitude2_Tensor Massive_Real_Subtraction::AddNLOTensor(METOOLS::Amplitude2_Tensor old_tensor){
  return old_tensor;
}


std::string Massive_Real_Subtraction::getType(){
  return "S";
}
