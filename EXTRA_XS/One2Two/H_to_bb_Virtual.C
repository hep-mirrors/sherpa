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
  SetVirtualMatrixFinite(momenta);
  SetVirtualMatrixE(momenta);
  SetVirtualMatrixE2(momenta);

  //Calculate(momenta, false);

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
          std::vector<std::pair<METOOLS::CSpinor<double>*, int>>> H_to_bb_Virtual::CalculateSpinors(const ATOOLS::Vec4D_Vector& momenta, bool anti){
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
      //for (DDSpin* sp : *v) {
      //  std::cout << "m_cur[1] hel="<<h<<" u0="<<(*sp)[0]<<" u1="<<(*sp)[1]<<" u2="<<(*sp)[2]<<" u3="<<(*sp)[3]<<"\n";
      //}
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
  M_finite = {
      M00, M01, M02, M03,
      M10, M11, M12, M13,
      M20, M21, M22, M23,
      M30, M31, M32, M33
    };
}


void H_to_bb_Virtual::SetVirtualMatrixE(const ATOOLS::Vec4D_Vector& momenta){
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
  M_epsilon = {
      M00, M01, M02, M03,
      M10, M11, M12, M13,
      M20, M21, M22, M23,
      M30, M31, M32, M33
    };
}


void H_to_bb_Virtual::SetVirtualMatrixE2(const ATOOLS::Vec4D_Vector& momenta){
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
  M_epsilon2 = {
      M00, M01, M02, M03,
      M10, M11, M12, M13,
      M20, M21, M22, M23,
      M30, M31, M32, M33
    };
}


std::map<std::string, std::complex<double>> H_to_bb_Virtual::CalculateBorn(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  typedef METOOLS::CSpinor<double> DDSpin;
  typedef std::pair<DDSpin*, int> SpinorWithHel;
  typedef std::vector<SpinorWithHel> SpinorVecWithHel;
  typedef std::pair<SpinorVecWithHel, SpinorVecWithHel> SpinorPairWithHel;
  using C = std::complex<double>;

  SpinorPairWithHel pair_spinors = CalculateSpinors(momenta, anti);

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

  std::map<std::string, C> born_res;
  born_res["00"]  = born_00;
  born_res["01"]  = born_01;
  born_res["10"]  = born_10;
  born_res["11"]  = born_11;

  return born_res;
}


std::map<std::string, std::map<std::string, std::complex<double>>> H_to_bb_Virtual::CalculateV(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
  // This method calculates the virtual amplitude for the H -> bb decay for the four helicity configurations
  // using the precomputed matrices and the spinors calculated from the currents.
  // There are finite terms, terms proportional to 1/epsilon and terms proportional to 1/epsilon^2.
  typedef METOOLS::CSpinor<double> DDSpin;
  typedef std::pair<DDSpin*, int> SpinorWithHel;
  typedef std::vector<SpinorWithHel> SpinorVecWithHel;
  typedef std::pair<SpinorVecWithHel, SpinorVecWithHel> SpinorPairWithHel;
  using C = std::complex<double>;

  SpinorPairWithHel pair_spinors = CalculateSpinors(momenta, anti);

  const SpinorVecWithHel &bottom = pair_spinors.first;
  const DDSpin* bottom_spinor_hel0 = (bottom[0].first); // first helicity state
  const DDSpin* bottom_spinor_hel1 = (bottom[1].first); // second helicity state
  const SpinorVecWithHel &antibottom = pair_spinors.second;
  const DDSpin* antibottom_spinor_hel0 = (antibottom[0].first); // first helicity state
  const DDSpin* antibottom_spinor_hel1 = (antibottom[1].first); // second helicity state

  if (!bottom_spinor_hel0 || !bottom_spinor_hel1 || !antibottom_spinor_hel0 || !antibottom_spinor_hel1) {
    THROW(fatal_error, "H_to_bb_Virtual::CalculateV - missing spinor(s) from currents");
  }

  // first matrix multiplication: M_finite * v(p_bbar) = "MV_ax" with: a = f, e or e2 for finite, 1/epsilon or 1/epsilon^2 and x = 0, 1 for helicity state
  std::vector<C> MV_f0(4, C(0.0,0.0));
  std::vector<C> MV_f1(4, C(0.0,0.0));
  std::vector<C> MV_e0(4, C(0.0,0.0));
  std::vector<C> MV_e1(4, C(0.0,0.0));
  std::vector<C> MV_e20(4, C(0.0,0.0));
  std::vector<C> MV_e21(4, C(0.0,0.0));

  // first matrix multiplication: M_finite * v(p_bbar)
  for (int i = 0; i <= 3; ++i) {
    for (int j = 0; j <= 3; ++j) {
      // M_finite is stored flat as 4x4 row-major: element (i,j) at index i*4 + j
      C M_f_ij = M_finite[i*4 + j];
      MV_f0[i] += M_f_ij * ((*antibottom_spinor_hel0)[j]);
      // 1/epsilon part
      C M_e_ij = M_epsilon[i*4 + j];
      MV_e0[i] += M_e_ij * ((*antibottom_spinor_hel0)[j]);
      // 1/epsilon^2 part
      C M_e2_ij = M_epsilon2[i*4 + j];
      MV_e20[i] += M_e2_ij * ((*antibottom_spinor_hel0)[j]);
    }
  }
  for (int i = 0; i <= 3; ++i) {
    for (int j = 0; j <= 3; ++j) {
      // M_finite is stored flat as 4x4 row-major: element (i,j) at index i*4 + j
      C M_f_ij = M_finite[i*4 + j];
      MV_f1[i] += M_f_ij * ((*antibottom_spinor_hel1)[j]);
      // 1/epsilon part
      C M_e_ij = M_epsilon[i*4 + j];
      MV_e1[i] += M_e_ij * ((*antibottom_spinor_hel1)[j]);
      // 1/epsilon^2 part
      C M_e2_ij = M_epsilon2[i*4 + j];
      MV_e21[i] += M_e2_ij * ((*antibottom_spinor_hel1)[j]);
    }
  }

  // second matrix multiplication: u(p_b) * (M_finite * v(p_bbar))
  // helicity configuration: 0,0
  C v_res_f00 = 0;
  C v_res_e00 = 0;
  C v_res_e200 = 0;
  for (int i = 0; i <= 3; ++i) {
    v_res_f00 += ((*bottom_spinor_hel0)[i]) * MV_f0[i];
    v_res_e00 += ((*bottom_spinor_hel0)[i]) * MV_e0[i];
    v_res_e200 += ((*bottom_spinor_hel0)[i]) * MV_e20[i];
  }
  // helicity configuration: 1,0
  C v_res_f10 = 0;
  C v_res_e10 = 0;
  C v_res_e210 = 0;
  for (int i = 0; i <= 3; ++i) {
    v_res_f10 += ((*bottom_spinor_hel1)[i]) * MV_f0[i];
    v_res_e10 += ((*bottom_spinor_hel1)[i]) * MV_e0[i];
    v_res_e210 += ((*bottom_spinor_hel1)[i]) * MV_e20[i];
  }
  // helicity configuration: 0,1
  C v_res_f01 = 0;
  C v_res_e01 = 0;
  C v_res_e201 = 0;
  for (int i = 0; i <= 3; ++i) {
    v_res_f01 += ((*bottom_spinor_hel0)[i]) * MV_f1[i];
    v_res_e01 += ((*bottom_spinor_hel0)[i]) * MV_e1[i];
    v_res_e201 += ((*bottom_spinor_hel0)[i]) * MV_e21[i];
  }
  // helicity configuration: 1,1
  C v_res_f11 = 0;
  C v_res_e11 = 0;
  C v_res_e211 = 0;
  for (int i = 0; i <= 3; ++i) {
    v_res_f11 += ((*bottom_spinor_hel1)[i]) * MV_f1[i];
    v_res_e11 += ((*bottom_spinor_hel1)[i]) * MV_e1[i];
    v_res_e211 += ((*bottom_spinor_hel1)[i]) * MV_e21[i];
  }

  // multiply with prefactor
  v_res_f00 *= VirtualPrefactor;
  v_res_e00 *= VirtualPrefactor;
  v_res_e200 *= VirtualPrefactor;

  v_res_f10 *= VirtualPrefactor;
  v_res_e10 *= VirtualPrefactor;
  v_res_e210 *= VirtualPrefactor;

  v_res_f01 *= VirtualPrefactor;
  v_res_e01 *= VirtualPrefactor;
  v_res_e201 *= VirtualPrefactor;

  v_res_f11 *= VirtualPrefactor;
  v_res_e11 *= VirtualPrefactor;
  v_res_e211 *= VirtualPrefactor;

  std::map<std::string, C> v_finite;
  v_finite["00"]  = v_res_f00;
  v_finite["01"]  = v_res_f01;
  v_finite["10"]  = v_res_f10;
  v_finite["11"]  = v_res_f11;

  std::map<std::string, C> v_epsilon;
  v_epsilon["00"]  = v_res_e00;
  v_epsilon["01"]  = v_res_e01;
  v_epsilon["10"]  = v_res_e10;
  v_epsilon["11"]  = v_res_e11;

  std::map<std::string, C> v_epsilon2;
  v_epsilon2["00"]  = v_res_e200;
  v_epsilon2["01"]  = v_res_e201;
  v_epsilon2["10"]  = v_res_e210;
  v_epsilon2["11"]  = v_res_e211;

  std::map<std::string, std::map<std::string, std::complex<double>>> res;
  res["v_finite"]  = std::move(v_finite);
  res["v_epsilon"] = std::move(v_epsilon);
  res["v_epsilon2"]= std::move(v_epsilon2);
  return res;
}


double H_to_bb_Virtual::CalculateVirtualCorrection(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  // Calculates the total virtual contribution with substraction terms already applied. Checks wether the epsilon terms cancel
  std::map<std::string, std::map<std::string, std::complex<double>>> virtual_amplitudes = CalculateV(momenta, anti);
  std::map<std::string, std::complex<double>> &v_finite  = virtual_amplitudes["v_finite"];
  std::map<std::string, std::complex<double>> &v_epsilon = virtual_amplitudes["v_epsilon"];
  std::map<std::string, std::complex<double>> &v_epsilon2= virtual_amplitudes["v_epsilon2"];

  std::map<std::string, std::complex<double>> born = CalculateBorn(momenta, anti); 

  // finite part
  std::complex<double> BV_f = born["00"] * std::conj(v_finite["00"]) + born["01"] * std::conj(v_finite["01"]) + born["10"] * std::conj(v_finite["10"]) + born["11"] * std::conj(v_finite["11"]);
  double v_correction_f = 2 * std::real(BV_f);

  // 1/epsilon part
  std::complex<double> BV_e = born["00"] * std::conj(v_epsilon["00"]) + born["01"] * std::conj(v_epsilon["01"]) + born["10"] * std::conj(v_epsilon["10"]) + born["11"] * std::conj(v_epsilon["11"]);
  double v_correction_e = 2 * std::real(BV_e);

  // 1/epsilon^2 part = 0 for H -> bb, because b is massive

  // substract (here: sum up) the substraction terms:
  double finite_sub = CalculateFiniteSubtraction(momenta);
  double epsilon_sub = CalculateEpsilonSubtraction(momenta);

  // todo: make sure that Epsilon terms cancel

  double v_correction = v_correction_f + finite_sub;
  double epsilon_sum = v_correction_e + epsilon_sub; // should be = 0"!!

  // print for test purpose
  std::cout << "finite_sub:" << finite_sub << std::endl;
  std::cout << "epsilon_sub:" << epsilon_sub << std::endl;
  std::cout << "v_correction_f:" << v_correction_f << std::endl;
  std::cout << "v_correction_e:" << v_correction_e << std::endl;
  std::cout << "Subtracted final results: " << std::endl;
  std::cout << "v_correction:" << v_correction << std::endl;
  std::cout << "v_correction_e:" << epsilon_sum << std::endl;

  return v_correction;
}


void H_to_bb_Virtual::SetUpPrefactors(const vector<Flavour>& flavs) {
  /* This method collects all constants that appear in the calculation and multiplies them.*/
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif
  double g_s = std::sqrt(4 * pi * alpha_qcd);
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

  double m_h = 125.09;
  double born_analytic_calc = 6 * (m_b/vev) * (m_b/vev) * (m_h*m_h - 4 * m_b*m_b);
  std::cout << "born_analytic_calc: " << born_analytic_calc << std::endl;
}


static double A(){
  // 1/epsilon prefactor in Gamma_j
  return 4.0/ 3.0;
}


static double B(double m_q, double mu){
  // epsilon-independent part in Gamma_j
  // input: quark mass and energy scale
  double C_F = 4.0/ 3.0;
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


static double nu_ijk_tilde(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  double m_j = std::sqrt(p_j * p_j);
  double m_k = std::sqrt(p_k * p_k);
  // either i or j is the gluon and massless. Therefore, the mass m_ij equals the bottom mass. Check, which mass is the bottom mass:
  double m_ij = m_j; // m_i = 0 (gluon), the emitter is only the bottom/ antibottom

  double mu_ij = mu_n(p_j, p_k, m_ij);
  double mu_k = mu_n(p_j, p_k, m_k);
  double lambda_val = lambda(1, mu_ij*mu_ij, mu_k*mu_k);
  return std::sqrt(lambda_val)/ (1 - mu_ij*mu_ij - mu_k*mu_k);
}


static double rho(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  double nu_tilde = nu_ijk_tilde(p_j, p_k);
  return std::sqrt((1 - nu_tilde) / (1 + nu_tilde));
}


static double C_j(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  // 1/epsilon prefactor in Nu_j
  return 1 / nu_jk(p_j, p_k) * std::log(rho(p_j, p_k));
}


static double Q_jk(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  double s_jk = 2 * p_j * p_k;
  double m2_j = p_j * p_j; // squared j mass
  double m2_k = p_k * p_k; // squared k mass
  double Q2_jk = s_jk + m2_j + m2_k; // squared Q_jk
  return std::sqrt(Q2_jk);
}


static double D_j(ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  // prefactor of 1/T_q^2 in Nu_j
  double C_F = 4.0 / 3.0;
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
  double var_nu_ijk_tilde = nu_ijk_tilde(p_j, p_k);
  double numerator = 1 - var_nu_ijk_tilde + 2*var_mu_n*var_mu_n / (1 - var_mu_j*var_mu_j - var_mu_k*var_mu_k);
  double denominator = 1 + var_nu_ijk_tilde + 2*var_mu_n*var_mu_n / (1 - var_mu_j*var_mu_j - var_mu_k*var_mu_k);
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
  double line2_2 = -2.0 * ( (var_Q_jk - m_k)*(var_Q_jk - m_k) - m_j) / (var_Q_jk*var_Q_jk);
  double line2_3 = 2.0*m_j*m_j / s_jk * std::log(m_j / (var_Q_jk - m_k));
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

  double C_F = 4.0 / 3.0;
  double gamma_j = 3.0 / 2.0 * C_F; // = gmma_q
  double s_jk = 2 * p_j * p_k;
  double K_j = (7.0/2.0 - pi*pi / 6)*C_F; // = K_q

  return gamma_j * std::log(mu*mu / s_jk) + gamma_j + K_j;
}


double H_to_bb_Virtual::CalculateFiniteSubtraction(const ATOOLS::Vec4D_Vector& momenta){
  // finite virtual substraction term
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif

  ATOOLS::Vec4<double> p_b = momenta[1];
  ATOOLS::Vec4<double> p_bb = momenta[2];

  double born_ME2 = 9.1018603234124864; // value taken from the H -> bbb calculation in Comix1to2, value right out of Decay_Channel::ME2(...)
                                        // => not multiplied with any symmetry factors/ colour factors
  double mu = 125.09;
  double alpha_qcd = MODEL::s_model -> ScalarFunction("alpha_S", mu*mu); // at Higgs scale
  double C_F = 4.0/3.0;
  double s_jk = 2 * p_b * p_bb;
  double gamma_E = 0.57721566490153286060; // Euler-Mascheroni constant
  double var_D_q = D_j(p_b, p_bb);
  double var_D_qq = D_j(p_bb, p_b);
  double var_E_q = E_j(p_b, p_bb);
  double var_E_qq = E_j(p_bb, p_b);
  double m_b = std::sqrt(p_b * p_b);

  double prefactor = born_ME2 * alpha_qcd / (2 * pi);
  double sum1 = (std::log(4 * pi * mu*mu / s_jk) - gamma_E) * (C_F * (C_j(p_b, p_bb) + C_j(p_bb, p_b)));
  double sum2 = var_D_q + var_D_qq + C_F * (var_E_q + var_E_qq - 2 * pi*pi*pi / 3);
  double sum3 = 2 * (A() * (std::log(4 * pi) - gamma_E) + B(m_b, mu) + F(p_b, p_bb, mu));

  return prefactor * (sum1 + sum2 + sum3);
}

double H_to_bb_Virtual::CalculateEpsilonSubtraction(const ATOOLS::Vec4D_Vector& momenta){
  // 1/epsilon prefactor of virtual substraction term
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif

  ATOOLS::Vec4<double> p_b = momenta[1];
  ATOOLS::Vec4<double> p_bb = momenta[2];

  double born_ME2 = 9.1018603234124864; // value taken from the H -> bbb calculation in Comix1to2, value right out of Decay_Channel::ME2(...)
                                        // => not multiplied with any symmetry factors/ colour factors
  double alpha_qcd = MODEL::s_model -> ScalarFunction("alpha_S", 125.09*125.09); // at Higgs scale
  double C_F = 4.0/3.0;

  return born_ME2 * alpha_qcd / (2 * pi) * (C_F * (C_j(p_b, p_bb) + C_j(p_bb, p_b)) + 2 * A());
}