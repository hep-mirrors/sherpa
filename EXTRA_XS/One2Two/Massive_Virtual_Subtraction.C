#include "EXTRA_XS/One2Two/Massive_Virtual_Subtraction.H"

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


Massive_Virtual_Subtraction::Massive_Virtual_Subtraction(const vector<Flavour>& flavs):
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(3), m_anticur(3), m_nhel(3), BornPrefactor(1.0)
{
  CalculateAlphaQCD(flavs[0].Mass()); // scale = mass of incoming particle

  if (flavs.size()!=3) THROW(fatal_error,"Internal error.");

  SetUpBornCurrents(flavs);
  SetUpBornPrefactor(flavs);
}


Massive_Virtual_Subtraction::~Massive_Virtual_Subtraction()
{
  for (size_t i(0);i<3;++i) {
    delete m_cur[i];
    delete m_anticur[i];
  }
  if (p_ci) delete p_ci;
}


void Massive_Virtual_Subtraction::CalculateAlphaQCD(double scale) {
  alpha_qcd = (MODEL::s_model) -> ScalarFunction("alpha_S", scale*scale); // at Higgs scale
}


size_t Massive_Virtual_Subtraction::NHel(const Flavour& fl) const
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


std::string Massive_Virtual_Subtraction::getType(){
  return "I";
}


void Massive_Virtual_Subtraction::SetUpBornPrefactor(const vector<Flavour>& flavs) {
  /* This method collects all constants (despite colour factors) that appear in the calculation and multiplies them.*/
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif
  double G_F = 1.16637886e-5; // Fermi constant in GeV^-2; value from PDG
  double vev = 1 / std::sqrt(G_F * std::sqrt(2)); // vacuum expectation value in GeV; doublecheck that value
  double m_b = flavs[2].Mass(); // b quark mass in GeV

  BornPrefactor = (-1) * m_b / vev;
}


void Massive_Virtual_Subtraction::SetUpBornCurrents(const vector<Flavour>& flavs){
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


void Massive_Virtual_Subtraction::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  CalculateBorn(momenta, anti); 
  double ME2_Born = 3 * std::real(born_hel["00"] * std::conj(born_hel["00"]) + born_hel["01"] * std::conj(born_hel["01"]) + born_hel["10"] * std::conj(born_hel["10"]) + born_hel["11"] * std::conj(born_hel["11"]));

  finite_sub = CalculateFiniteSubtraction(momenta, ME2_Born);
  double epsilon_sub = CalculateEpsilonSubtraction(momenta, ME2_Born);
  // todo: fill this spin amplitue object with helicity-dependent values of I
}


double Massive_Virtual_Subtraction::get_NLO_ME2(){
  return finite_sub;
}


std::pair<std::vector<std::pair<METOOLS::CSpinor<double>*, int>>,
          std::vector<std::pair<METOOLS::CSpinor<double>*, int>>> Massive_Virtual_Subtraction::CalculateBornSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti){
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


void Massive_Virtual_Subtraction::CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti){
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
    THROW(fatal_error, "Massive_Virtual_Subtraction::CalculateV - missing spinor(s) from currents");
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
}


static double A(){
  // 1/epsilon prefactor in Gamma_j; C_F
  double C_F = 4.0/3.0;
  return C_F;
}


static double B(double m_q, double mu){
  // epsilon-independent part in Gamma_j
  // input: quark mass and energy scale
  double C_F = 4.0/3.0;
  return C_F * (1/2 * std::log(m_q*m_q/ mu*mu) - 2); 
}


static double nu_jk(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  // relative velocity between two massive momenta p_j and p_k
  return std::sqrt(1 - (p_j*p_j)*(p_k*p_k) / ((p_j*p_k) * (p_j*p_k)));
}


static double mu_n(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k, double m_n){
  // rescaled parton masses
  ATOOLS::Vec4<double> Q = p_j+ p_k; // total outgiong momentum
  return m_n / (std::sqrt(Q*Q));
}


static double lambda(double x, double y, double z){
  return x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
}


static double nu_jk_tilde(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  double m_j = std::sqrt(p_j * p_j);
  double m_k = std::sqrt(p_k * p_k);

  double mu_j = mu_n(p_j, p_k, m_j);
  double mu_k = mu_n(p_j, p_k, m_k);
  double lambda_val = lambda(1, mu_j*mu_j, mu_k*mu_k);
  return std::sqrt(lambda_val)/ (1 - mu_j*mu_j - mu_k*mu_k);
}


static double rho(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  double nu_tilde = nu_jk_tilde(p_j, p_k);
  return std::sqrt((1 - nu_tilde) / (1 + nu_tilde));
}


static double C_j(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  // 1/epsilon prefactor in Nu_j
  return 1 / nu_jk(p_j, p_k) * std::log(rho(p_j, p_k));
}


static double Q_jk(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  ATOOLS::Vec4<double> Q_jk = p_j + p_k;
  double Q2_jk = (p_j + p_k) * (p_j + p_k);

  return std::sqrt(Q2_jk);
}


static double D_j(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  // prefactor of 1/T_q^2 in Nu_j
  double C_F = 4.0/3.0;
  double s_jk = 2 * p_j * p_k;
  double gamma_q = 3.0 / 2.0 * C_F;
  return gamma_q * std::log(s_jk/ (Q_jk(p_j, p_k) * Q_jk(p_j, p_k)));
}


static double rho_n(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k, ATOOLS::Vec4<double> p_n){
  double m_n = std::sqrt(p_n*p_n);
  double m_j = std::sqrt(p_j*p_j);
  double m_k = std::sqrt(p_k*p_k);

  double var_mu_n = mu_n(p_j, p_k, m_n);
  double var_mu_j = mu_n(p_j, p_k, m_j);
  double var_mu_k = mu_n(p_j, p_k, m_k);
  double var_nu_jk_tilde = nu_jk_tilde(p_j, p_k);
  double numerator = 1 - var_nu_jk_tilde + 2*var_mu_n*var_mu_n / (1 - var_mu_j*var_mu_j - var_mu_k*var_mu_k);
  double denominator = 1 + var_nu_jk_tilde + 2*var_mu_n*var_mu_n / (1 - var_mu_j*var_mu_j - var_mu_k*var_mu_k);
  return std::sqrt(numerator / denominator);
}


static double li2(double x){
  // function taken from dilog.C

  // routines only valid for real values of the argument; 
  // imaginary parts not presently given.
  // this version uses 't Hooft and Veltman's change of variable
  // good to ~ 10^(-16)
  const double PISQ6  =  1.64493406684822643647;
  double x_0 = -0.30;
  double x_1 = 0.25;
  double x_2 = 0.51;
  if (x == 1.) return PISQ6;
  if (x <= x_0){ 
    double temp = std::log(Abs(1.0-x));
    return -li2(-x/(1.0-x)) - temp*temp/2 ; }
  else if (x < x_1){
    double z = - std::log(1.0-x);
    double temp = z*(1.0-z/4.0*(1.0-z/9.0*(1.0-z*z/100.0
                  *(1.0-5.0*z*z/294.0*(1.0-7.0*z*z/360.0
                  *(1.0-5.0*z*z/242.0*(1.0-7601.0*z*z/354900.0
                  *(1.0-91.0*z*z/4146.0*(1.0-3617.0*z*z/161840.0)
                   ))))))));
    return temp; }
    else if (x < x_2) return - li2(-x) + li2(x*x)/2.0 ;
    else { return PISQ6 - li2(1.0-x) 
                  - std::log(Abs(x))*std::log(Abs(1.0-x)) ; }
}


static double E_j(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  // collects the other terms in Nu_j

  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif

  // pre-calculate some variables
  double var_rho = rho(p_j, p_k);
  double var_rho_j = rho_n(p_j, p_k, p_j);
  double var_rho_k = rho_n(p_j, p_k, p_k);
  double var_Q_jk = Q_jk(p_j, p_k);
  double s_jk = 2 * p_j * p_k;
  double m_j = std::sqrt(p_j * p_j);
  double m_k = std::sqrt(p_k * p_k);

  // Nu^s - part (without the 1/epsilon term):
  double sum1 = -1.0/4.0 * std::log(var_rho_j*var_rho_j) * std::log(var_rho_j*var_rho_j);
  double sum2 = -1.0/4.0 * std::log(var_rho_k*var_rho_k) * std::log(var_rho_k*var_rho_k);
  double sum3 = 1/nu_jk(p_j, p_k) * std::log(var_rho)* std::log(var_Q_jk*var_Q_jk/ s_jk);
  double part1 = 1/nu_jk(p_j, p_k) * (sum1 + sum2 - pi*pi/6.0) + sum3;

  // Nu^NS - part (without the 1/T_q^2 term):
  // this expression is very long and therefore sorted according to the lines in the Catani Dittmaier paper (formula 6.21)
  double line1_1 = std::log(var_rho * var_rho) * std::log(1 + var_rho * var_rho);
  double line1_2 = 2 * li2(var_rho * var_rho);
  double line1_3 = - li2(1 - var_rho_j * var_rho_j);
  double line1_4 = - li2(1 - var_rho_k * var_rho_k);
  double line1 = 1/nu_jk(p_j, p_k) * (line1_1 + line1_2 + line1_3 + line1_4 - (pi*pi)/6.0);  

  double line2_1 = std::log((var_Q_jk - m_k) / var_Q_jk);
  double line2_2 = -2.0 * std::log(((var_Q_jk - m_k)*(var_Q_jk - m_k) - m_j) / (var_Q_jk*var_Q_jk));
  double line2_3 = -2.0*m_j*m_j / s_jk * std::log(m_j / (var_Q_jk - m_k));
  double line2 = line2_1 + line2_2 + line2_3;

  double line3_1 = -m_k / (var_Q_jk - m_k);
  double line3_2 = (2.0 * m_k * (2 * m_k - var_Q_jk)) / s_jk;
  double line3_3 = pi*pi/ 2.0;
  double line3 = line3_1 + line3_2 + line3_3;
  return part1 + line1 + line2 + line3;
}


static double F(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k, double mu){
  // other terms collected
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif
  double C_F = 4.0/3.0;
  double gamma_j = 3.0 / 2.0 * C_F; // = gmma_q
  double s_jk = 2 * p_j * p_k;
  double K_q = (7.0/2.0 - pi*pi / 6)*C_F; // = K_q

  return gamma_j * std::log(mu*mu / s_jk) + gamma_j + K_q;
}


double Massive_Virtual_Subtraction::CalculateFiniteSubtraction(const ATOOLS::Vec4D_Vector& momenta, double born_ME2){
  // finite virtual subtraction term
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif

  ATOOLS::Vec4<double> p_b = momenta[1];
  ATOOLS::Vec4<double> p_bb = momenta[2];
  double C_F = 4.0/3.0;

  double mu = std::sqrt( (p_b + p_bb) * (p_b + p_bb));
  double s_jk = 2 * p_b * p_bb;
  double gamma_E = 0.57721566490153286060; // Euler-Mascheroni constant
  double var_D_q = D_j(p_b, p_bb);
  double var_D_qq = D_j(p_bb, p_b);
  double var_E_q = E_j(p_b, p_bb);
  double var_E_qq = E_j(p_bb, p_b);
  double m_b = std::sqrt(p_b * p_b);

  double prefactor = alpha_qcd / (2 * pi);
  double sum1 = (std::log(4 * pi * mu*mu / s_jk) - gamma_E) * (C_F * (C_j(p_b, p_bb) + C_j(p_bb, p_b)));
  double sum2 = var_D_q + var_D_qq + C_F * (var_E_q + var_E_qq - 2 * pi*pi*pi / 3);
  double sum3 = 2 * (A() * (std::log(4 * pi) - gamma_E) + B(m_b, mu) + F(p_b, p_bb, mu));

  for (size_t i=0; i<size(); ++i) { // reset values of spin_amplitude object
     (*this)[i] = Complex(0.0, 0.0);
  }
  
  (*this)[0] = born_hel["00"];
  (*this)[1] = born_hel["10"];
  (*this)[2] = born_hel["01"];
  (*this)[3] = born_hel["11"];

  for (size_t i=0; i<size(); ++i) {
   (*this)[i] *= std::sqrt(prefactor * (sum1 + sum2 + sum3));
  }

  return born_ME2 * prefactor * (sum1 + sum2 + sum3);
}


double Massive_Virtual_Subtraction::CalculateEpsilonSubtraction(const ATOOLS::Vec4D_Vector& momenta, double born_ME2){
  // 1/epsilon prefactor of virtual subtraction term
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif

  ATOOLS::Vec4<double> p_b = momenta[1];
  ATOOLS::Vec4<double> p_bb = momenta[2];
  double C_F = 4.0/3.0;

  return born_ME2 * alpha_qcd / (2 * pi) * (C_F * (C_j(p_b, p_bb) + C_j(p_bb, p_b)) + 2 * A());
}
