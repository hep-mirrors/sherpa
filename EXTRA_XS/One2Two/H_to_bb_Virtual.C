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


std::map<std::string, std::complex<double>> H_to_bb_Virtual::CalculateV(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
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

  std::map<std::string, C> vres;
  vres["v_res_f00"]  = v_res_f00;   vres["v_res_e00"]  = v_res_e00;   vres["v_res_e200"]  = v_res_e200;
  vres["v_res_f10"]  = v_res_f10;   vres["v_res_e10"]  = v_res_e10;   vres["v_res_e210"]  = v_res_e210;
  vres["v_res_f01"]  = v_res_f01;   vres["v_res_e01"]  = v_res_e01;   vres["v_res_e201"]  = v_res_e201;
  vres["v_res_f11"]  = v_res_f11;   vres["v_res_e11"]  = v_res_e11;   vres["v_res_e211"]  = v_res_e211;

  std::map<std::string, C> v_finite;
  v_finite["00"]  = v_res_f00;
  v_finite["10"]  = v_res_f10;
  v_finite["01"]  = v_res_f01;
  v_finite["11"]  = v_res_f11;

  return v_finite;
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