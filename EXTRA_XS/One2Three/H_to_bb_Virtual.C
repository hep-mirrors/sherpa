#include "EXTRA_XS/One2Three/H_to_bb_Virtual.H"

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

H_to_bb_Virtual::H_to_bb_Virtual(const vector<Flavour>& flavs, MODEL::Model_Base* s_model):
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(3), m_anticur(3), m_nhel(3), 
  BornPrefactor(1.0), VirtualPrefactor(1.0)
{
  Calculate_alpha_QCD(s_model);

  if (flavs.size()!=3) THROW(fatal_error,"Internal error.");

  // example momenta
  std::vector<Vec4D> momenta(3);
  momenta[0]=Vec4D(125.0,0.0,0.0,0.0); // Higgs
  momenta[1]=Vec4D(62.54499999999998, 2.3910658074969522, 43.056947854781868, 45.033905790357714); // b
  momenta[2]=Vec4D(62.545000000000009, -2.3910658074969575, -43.056947854781889, -45.033905790357707); // bbar

  SetUpCurrents(flavs);
  SetUpPrefactors(flavs);
  CalculateBorn(momenta);

  // test:
  std::array<double, 4> p1 = {0.0, 0.0, 0.0, 0.0};
  std::array<double, 4> p2 = {0.0, 0.0, 0.0, 0.0};
  std::array<std::complex<double>,16> M_finite(GetVirtualMatrixFinite(momenta));
  std::array<std::complex<double>,16> M_epsilon(GetVirtualMatrixE(momenta));
  std::array<std::complex<double>,16> M_epsilon2(GetVirtualMatrixE2(momenta));

  // Ausgabe aller Elemente
  std::cout << "\n=== M_finite Matrix ===" << std::endl;
  std::cout << M_finite[0] << " " << M_finite[1] << " " << M_finite[2] << " " << M_finite[3] << std::endl;
  std::cout << M_finite[4] << " " << M_finite[5] << " " << M_finite[6] << " " << M_finite[7] << std::endl;
  std::cout << M_finite[8] << " " << M_finite[9] << " " << M_finite[10] << " " << M_finite[11] << std::endl;
  std::cout << M_finite[12] << " " << M_finite[13] << " " << M_finite[14] << " " << M_finite[15] << std::endl;

}

H_to_bb_Virtual::~H_to_bb_Virtual()
{
  for (size_t i(0);i<3;++i) {
    delete m_cur[i];
    delete m_anticur[i];
  }
}

void H_to_bb_Virtual::Calculate_alpha_QCD(MODEL::Model_Base* s_model) {
  alpha_qcd = s_model -> ScalarFunction("alpha_S", 15625); // at Higgs scale
  std::cout << "The cpl value is: " << alpha_qcd << std::endl;
}


size_t H_to_bb_Virtual::NHel(const Flavour& fl)
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


void H_to_bb_Virtual::CalculateBorn(const ATOOLS::Vec4D_Vector& momenta){
  p_ci->GeneratePoint();
  // analytical try, just current, not anticurrent

  // copy currents
  std::vector<METOOLS::Current*> born_mcur(m_cur); 
  METOOLS::Current* bottom_cur = born_mcur[1];
  METOOLS::Current* antibottom_cur = born_mcur[1];

  // copy anticurrents
  std::vector<METOOLS::Current*> born_anticur(m_anticur); 
  METOOLS::Current* bottom_anticur = born_anticur[1];
  METOOLS::Current* antibottom_anticur = born_anticur[1];

  for (size_t i(0);i<m_cur.size();++i) {
    m_cur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci->I()[i],p_ci->J()[i],0);
    m_cur[i]->Print();
  }

  //vector<int> fill(m_n,1); // output amplitude vector
  //for (size_t i(0);i<m_n;++i) (*this)[i]=Complex(0.0,0.0);
  CF<double>* b_fermion = dynamic_cast<CF<double>*>(m_cur[1]);
  //CF<double>* bbar_fermion = dynamic_cast<CF<double>*>(m_cur[2]);

  typedef METOOLS::CSpinor<double> DDSpin;
  DDSpin *spin = DDSpin::New();  
  double m2 = 25.0;
  const double helicity = 1;
  const int mass_sign = 1;
  spin->Construct(helicity, momenta[1], m2, mass_sign);

  // Zugriff auf m_u über operator[]
  std::complex<double> u0 = (*spin)[0];
  std::complex<double> u1 = (*spin)[1];
  std::complex<double> u2 = (*spin)[2];
  std::complex<double> u3 = (*spin)[3];

  std::cout << "CSpinor components: "
            << u0 << " " << u1 << " " << u2 << " " << u3 << std::endl;
}


std::array<std::complex<double>,16> H_to_bb_Virtual::GetVirtualMatrixFinite(ATOOLS::Vec4D_Vector& momenta) {
    /* This method provides precomputed parts of the finite virtual correction matrix element for the H -> bb decay.
   * 
   * The matrix correspond to the expression: gamma-matrix * Loop-Integral * gamma-matrix.
   * 
   * Since the virtual corrections are momentum-independent for this process, the calculation
   * could be performed externally using a Python script.
   * 
   * The matrix contains three components:
   * - Finite part: The UV-finite contribution to the virtual amplitude
   * - 1/epsilon and 1/epsilon^2 divergence: The dimensional regularization poles that will be cancelled
   *   by substraction terms
   * 
   * Here is the original Python code located: 
   * https://github.com/LeaBaumann/pre-calculations_integrals_and_subtraction/blob/main/Integrals_and_gammas.py
   */
  
  using C = std::complex<double>;

  // Calculate q1 and q2: q1 = p_bbar; q2 = p_bbar + p_Higgs
    std::array<double, 4> q1 = {momenta[2][0], momenta[2][1], momenta[2][2], momenta[2][3]};
    
    std::array<double, 4> q2 = {
        momenta[2][0] + momenta[0][0],
        momenta[2][1] + momenta[0][1], 
        momenta[2][2] + momenta[0][2],
        momenta[2][3] + momenta[0][3]
    };

  const double diag_re = -0.002618352*q1[0]*q1[0] + 0.005900148*q1[0]*q2[0] + 0.002618352*q1[1]*q1[1] - 0.005900148*q1[1]*q2[1]
      + 0.002618352*q1[2]*q1[2] - 0.005900148*q1[2]*q2[2] + 0.002618352*q1[3]*q1[3] - 0.005900148*q1[3]*q2[3]
      + 0.000698944*q2[0]*q2[0] - 0.000698944*q2[1]*q2[1] - 0.000698944*q2[2]*q2[2] - 0.000698944*q2[3]*q2[3]
      - 26.322871021248;

  const double diag_im = + 0.001209628*q1[0]*q1[0] - 0.006166716504*q1[0]*q2[0] - 0.001209628*q1[1]*q1[1] + 0.006166716504*q1[1]*q2[1]
      - 0.001209628*q1[2]*q1[2] + 0.006166716504*q1[2]*q2[2] - 0.001209628*q1[3]*q1[3] + 0.006166716504*q1[3]*q2[3]
      - 0.000401548*q2[0]*q2[0] + 0.000401548*q2[1]*q2[1] + 0.000401548*q2[2]*q2[2] + 0.000401548*q2[3]*q2[3]
      + 12.33920955744;
  
  const C d(diag_re, diag_im);

  // elements of the matrix: 
  C M00 = d;
  C M01 = C(0.0, 0.0);
  C M02 = C(+0.00368220816*q1[0] - 0.00368220816*q1[3] - 0.01387288464*q2[0] + 0.01387288464*q2[3],
              -7.57929599999995e-5*q1[0] + 7.57929599999995e-5*q1[3] + 0.01516399104*q2[0] - 0.01516399104*q2[3]);
  C M03 = C(-0.00368220816*q1[1] + 7.57929599999995e-5*q1[2] + 0.01387288464*q2[1] - 0.01516399104*q2[2],
               7.57929599999995e-5*q1[1] + 0.00368220816*q1[2] - 0.01516399104*q2[1] - 0.01387288464*q2[2]);

  C M10 = C(0.0, 0.0);
  C M11 = d;
  C M12 = C(-0.00368220816*q1[1] - 7.57929599999995e-5*q1[2] + 0.01387288464*q2[1] + 0.01516399104*q2[2],
               7.57929599999995e-5*q1[1] - 0.00368220816*q1[2] - 0.01516399104*q2[1] + 0.01387288464*q2[2]);
  C M13 = C(+0.00368220816*q1[0] + 0.00368220816*q1[3] - 0.01387288464*q2[0] - 0.01387288464*q2[3],
              -7.57929599999995e-5*q1[0] - 7.57929599999995e-5*q1[3] + 0.01516399104*q2[0] + 0.01516399104*q2[3]);

  C M20 = C(+0.00368220816*q1[0] + 0.00368220816*q1[3] - 0.01387288464*q2[0] - 0.01387288464*q2[3],
              -7.57929599999995e-5*q1[0] - 7.57929599999995e-5*q1[3] + 0.01516399104*q2[0] + 0.01516399104*q2[3]);
  C M21 = C(+0.00368220816*q1[1] - 7.57929599999995e-5*q1[2] - 0.01387288464*q2[1] + 0.01516399104*q2[2],
              -7.57929599999995e-5*q1[1] - 0.00368220816*q1[2] + 0.01516399104*q2[1] + 0.01387288464*q2[2]);
  C M22 = d;
  C M23 = C(0.0, 0.0);

  C M30 = C(+0.00368220816*q1[1] + 7.57929599999995e-5*q1[2] - 0.01387288464*q2[1] - 0.01516399104*q2[2],
              -7.57929599999995e-5*q1[1] + 0.00368220816*q1[2] + 0.01516399104*q2[1] - 0.01387288464*q2[2]);
  C M31 = C(+0.00368220816*q1[0] - 0.00368220816*q1[3] - 0.01387288464*q2[0] + 0.01387288464*q2[3],
              -7.57929599999995e-5*q1[0] + 7.57929599999995e-5*q1[3] + 0.01516399104*q2[0] - 0.01516399104*q2[3]);
  C M32 = C(0.0, 0.0);
  C M33 = d;

  // put matrix together
  std::array<C,16> M = {
      M00, M01, M02, M03,
      M10, M11, M12, M13,
      M20, M21, M22, M23,
      M30, M31, M32, M33
    };

  return M;
}


std::array<std::complex<double>,16> H_to_bb_Virtual::GetVirtualMatrixE(ATOOLS::Vec4D_Vector& momenta){
    /* This method provides precomputed parts of the finite virtual correction matrix element for the H -> bb decay.
   * 
   * The matrix correspond to the expression: gamma-matrix * Loop-Integral * gamma-matrix.
   * 
   * Since the virtual corrections are momentum-independent for this process, the calculation
   * could be performed externally using a Python script.
   * 
   * The matrix contains three components:
   * - Finite part: The UV-finite contribution to the virtual amplitude
   * - 1/epsilon and 1/epsilon^2 divergence: The dimensional regularization poles that will be cancelled
   *   by substraction terms
   * 
   * Here is the original Python code located: 
   * https://github.com/LeaBaumann/pre-calculations_integrals_and_subtraction/blob/main/Integrals_and_gammas.py
   */
  using C = std::complex<double>;

  // Calculate q1 and q2: q1 = p_bbar; q2 = p_bbar + p_Higgs
    std::array<double, 4> q1 = {momenta[2][0], momenta[2][1], momenta[2][2], momenta[2][3]};
    
    std::array<double, 4> q2 = {
        momenta[2][0] + momenta[0][0],
        momenta[2][1] + momenta[0][1], 
        momenta[2][2] + momenta[0][2],
        momenta[2][3] + momenta[0][3]
    };

      // Calculate diagonal elements
  const double diag_re = -0.001658656*q1[0]*q2[0] + 0.001658656*q1[1]*q2[1] + 0.001658656*q1[2]*q2[2] + 0.001658656*q1[3]*q2[3] - 0.0401500905984;
  const double diag_im = 0.000805588*q1[0]*q2[0] - 0.000805588*q1[1]*q2[1] - 0.000805588*q1[2]*q2[2] - 0.000805588*q1[3]*q2[3] + 0.0195003853632;
  
  const C d(diag_re, diag_im);

  // elements of the matrix: 
  C M00 = d;
  C M01 = C(0.0, 0.0);
  C M02 = C(0.000829328*q1[0] - 0.000829328*q1[3] + 0.00408029376*q2[0] - 0.00408029376*q2[3],
            -0.000402794*q1[0] + 0.000402794*q1[3] - 0.00198174648*q2[0] + 0.00198174648*q2[3]);
  C M03 = C(-0.000829328*q1[1] + 0.000402794*q1[2] - 0.00408029376*q2[1] + 0.00198174648*q2[2],
            0.000402794*q1[1] + 0.000829328*q1[2] + 0.00198174648*q2[1] + 0.00408029376*q2[2]);

  C M10 = C(0.0, 0.0);
  C M11 = d;
  C M12 = C(-0.000829328*q1[1] - 0.000402794*q1[2] - 0.00408029376*q2[1] - 0.00198174648*q2[2],
            0.000402794*q1[1] - 0.000829328*q1[2] + 0.00198174648*q2[1] - 0.00408029376*q2[2]);
  C M13 = C(0.000829328*q1[0] + 0.000829328*q1[3] + 0.00408029376*q2[0] + 0.00408029376*q2[3],
            -0.000402794*q1[0] - 0.000402794*q1[3] - 0.00198174648*q2[0] - 0.00198174648*q2[3]);

  C M20 = C(0.000829328*q1[0] + 0.000829328*q1[3] + 0.00408029376*q2[0] + 0.00408029376*q2[3],
            -0.000402794*q1[0] - 0.000402794*q1[3] - 0.00198174648*q2[0] - 0.00198174648*q2[3]);
  C M21 = C(0.000829328*q1[1] - 0.000402794*q1[2] + 0.00408029376*q2[1] - 0.00198174648*q2[2],
            -0.000402794*q1[1] - 0.000829328*q1[2] - 0.00198174648*q2[1] - 0.00408029376*q2[2]);
  C M22 = d;
  C M23 = C(0.0, 0.0);

  C M30 = C(0.000829328*q1[1] + 0.000402794*q1[2] + 0.00408029376*q2[1] + 0.00198174648*q2[2],
            -0.000402794*q1[1] + 0.000829328*q1[2] - 0.00198174648*q2[1] + 0.00408029376*q2[2]);
  C M31 = C(0.000829328*q1[0] - 0.000829328*q1[3] + 0.00408029376*q2[0] - 0.00408029376*q2[3],
            -0.000402794*q1[0] + 0.000402794*q1[3] - 0.00198174648*q2[0] + 0.00198174648*q2[3]);
  C M32 = C(0.0, 0.0);
  C M33 = d;

  // put matrix together
  std::array<C,16> M = {
      M00, M01, M02, M03,
      M10, M11, M12, M13,
      M20, M21, M22, M23,
      M30, M31, M32, M33
    };

  return M;
}


std::array<std::complex<double>,16> H_to_bb_Virtual::GetVirtualMatrixE2(ATOOLS::Vec4D_Vector& momenta){
    /* This method provides precomputed parts of the finite virtual correction matrix element for the H -> bb decay.
   * 
   * The matrix correspond to the expression: gamma-matrix * Loop-Integral * gamma-matrix.
   * 
   * Since the virtual corrections are momentum-independent for this process, the calculation
   * could be performed externally using a Python script.
   * 
   * The matrix contains three components:
   * - Finite part: The UV-finite contribution to the virtual amplitude
   * - 1/epsilon and 1/epsilon^2 divergence: The dimensional regularization poles that will be cancelled
   *   by substraction terms
   * 
   * Here is the original Python code located: 
   * https://github.com/LeaBaumann/pre-calculations_integrals_and_subtraction/blob/main/Integrals_and_gammas.py
   */
  
  using C = std::complex<double>;

  // Calculate q1 and q2: q1 = p_bbar; q2 = p_bbar + p_Higgs
    std::array<double, 4> q1 = {momenta[2][0], momenta[2][1], momenta[2][2], momenta[2][3]};
    
    std::array<double, 4> q2 = {
        momenta[2][0] + momenta[0][0],
        momenta[2][1] + momenta[0][1], 
        momenta[2][2] + momenta[0][2],
        momenta[2][3] + momenta[0][3]
    };

    // elements of the matrix: 
  C M00 = C(0.0, 0.0);
  C M01 = C(0.0, 0.0);
  C M02 = C(0.0, 0.0);
  C M03 = C(0.0, 0.0);

  C M10 = C(0.0, 0.0);
  C M11 = C(0.0, 0.0);
  C M12 = C(0.0, 0.0);
  C M13 = C(0.0, 0.0);

  C M20 = C(0.0, 0.0);
  C M21 = C(0.0, 0.0);
  C M22 = C(0.0, 0.0);
  C M23 = C(0.0, 0.0);

  C M30 = C(0.0, 0.0);
  C M31 = C(0.0, 0.0);
  C M32 = C(0.0, 0.0);
  C M33 = C(0.0, 0.0);

  // put matrix together
  std::array<C,16> M = {
      M00, M01, M02, M03,
      M10, M11, M12, M13,
      M20, M21, M22, M23,
      M30, M31, M32, M33
    };
  
  return M;
}


void H_to_bb_Virtual::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
} 


void H_to_bb_Virtual::SetUpPrefactors(const vector<Flavour>& flavs) {
  /* This method collects all constants that appear in the calculation and multiplies them.
  There are constante both in the real and the virtual correction. */
  double g_s = std::sqrt(4 * std::acos(-1) * alpha_qcd);  // strong gauge coupling; std::acos(-1) = pi
  double G_F = 1.16637886e-5; // Fermi constant in GeV^-2; value from PDG
  double vev = 1 / std::sqrt(G_F * std::sqrt(2)); // vacuum expectation value in GeV; doublecheck that value
  double m_b = flavs[2].Mass(); // b quark mass in GeV
  std::cout << "The vev value is: " << vev << std::endl;

  //double mu = flavs[0].Mass(); // Renormalisation scale: Higgs mass in GeV
  //double epsilon = 0.0001; // small parameter for dimensional regularization

  VirtualPrefactor = (-1) * std::pow(g_s, 2) * m_b / vev;
  std::cout << "Virtual prefactor: " << VirtualPrefactor << std::endl;

  BornPrefactor = (-1) * m_b / vev;
  std::cout << "Born prefactor: " << BornPrefactor << std::endl;
}