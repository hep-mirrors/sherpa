#include "EXTRA_XS/One2Three/H_to_bb_Virtual.H"

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "EXTRA_XS/Main/ME_Tools.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include <assert.h>

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace PHASIC;
using namespace std;
using namespace MODEL;

H_to_bb_Virtual::H_to_bb_Virtual(const vector<Flavour>& flavs, const Flavour& prop1, const Flavour& prop2,
                    MODEL::Model_Base* s_model) :
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(2), m_anticur(2), m_nhel(2),
  m_prop(prop1) // change this
{
  Calculate_alpha_QCD(s_model);

  //DEBUG_FUNC(flavs<<" with prop "<<prop<<" in "<<gluon<<","<<propj);
  //assert(non_prop>0 && gluon>0 && propj>0);
  //if (flavs.size()!=4) THROW(fatal_error,"Internal error."); // change later
  Vec4D k(1.0,0.0,1.0,0.0); // gauge 

  // flavs[1]: gluon; flavs[2]: b quark; flavs[3]: bbar quark
  Current_Key p1_cur(flavs[2],MODEL::s_model,1);
  m_cur[0] = Current_Getter::GetObject("D"+p1_cur.Type(),p1_cur); // get the object that corresponds to the key "u_bar"
  if (m_cur[0]==NULL) THROW(fatal_error, "current not found");
  m_cur[0]->SetDirection(-1); // +1 for i=0, -1 for i=1, 2, 3
  m_cur[0]->SetId(std::vector<int>(1,0));
  m_cur[0]->InitPols(std::vector<int>(1,m_spins[0])); // define allowed polarization states
  m_cur[0]->SetKey(0);
  m_cur[0]->SetGauge(k);
  m_nhel[0]=NHel(flavs[2]); // number of helicity states (based on spin properties)

  Current_Key p2_cur(flavs[3],MODEL::s_model,1);
  m_cur[1] = Current_Getter::GetObject("D"+p2_cur.Type(),p2_cur); // get the object that corresponds to the key "ckey"
  if (m_cur[1]==NULL) THROW(fatal_error, "current not found");
  m_cur[1]->SetDirection(-1); // +1 for i=0, -1 for i=1, 2, 3
  m_cur[1]->SetId(std::vector<int>(1,1));
  m_cur[1]->InitPols(std::vector<int>(1,m_spins[1])); // define allowed polarization states
  m_cur[1]->SetKey(1);
  m_cur[1]->SetGauge(k);
  m_nhel[3]=NHel(flavs[1]); // number of helicity states (based on spin properties)

// do the same for the anticurrent
  Current_Key p1_anticur(flavs[2].Bar(),MODEL::s_model,1);
  m_anticur[0] = Current_Getter::GetObject("D"+p1_anticur.Type(),p1_anticur);
  if (m_anticur[0]==NULL) THROW(fatal_error, "current not found");
  m_anticur[0]->SetDirection(-1);
  m_anticur[0]->SetId(std::vector<int>(1,0));
  m_anticur[0]->InitPols(std::vector<int>(1,m_spins[0]));
  m_anticur[0]->SetKey(0);
  m_anticur[0]->SetGauge(k);

  Current_Key p2_anticur(flavs[3].Bar(),MODEL::s_model,1);
  m_anticur[1] = Current_Getter::GetObject("D"+p2_anticur.Type(),p2_anticur);
  if (m_anticur[1]==NULL) THROW(fatal_error, "current not found");
  m_anticur[1]->SetDirection(-1);
  m_anticur[1]->SetId(std::vector<int>(1,1));
  m_anticur[1]->InitPols(std::vector<int>(1,m_spins[1]));
  m_anticur[1]->SetKey(1);
  m_anticur[1]->SetGauge(k);

  Calculate_ME2(flavs);
}

H_to_bb_Virtual::~H_to_bb_Virtual()
{
  for (size_t i(0);i<4;++i) {
    delete m_cur[i];
    delete m_anticur[i];
  }
  delete m_scur;
  delete m_antiscur;
  delete m_fcur;
  delete m_antifcur;
  if (p_ci) delete p_ci;
}

void H_to_bb_Virtual::Calculate_alpha_QCD(MODEL::Model_Base* s_model) {
  /*s_model->GetCouplings(m_cpls);
  if (m_cpls.find("Alpha_QCD") != m_cpls.end())
    p_aqcd=m_cpls.Get("Alpha_QCD");
  if (p_aqcd)
    std::cout << "p_aqcd: " << *p_aqcd << std::endl;
  else
    std::cout << "p_aqcd is null." << std::endl;

  // this does not really work yet (setting new scale, because here we have alpha_qcd scaled with Z mass):
  double Scale = p_aqcd->Scale();
  std::cout << "Scale: " << Scale << std::endl;
  double newScale = 15625; // Higgs squared
  p_aqcd->SetScale(&newScale);
  p_aqcd->Calculate();
  std::cout << "New Scale: " << p_aqcd->Scale() << std::endl;

  alpha_qcd = p_aqcd->Default(); // at Z scale
  std::cout << "The cpl value is: " << alpha_qcd << std::endl;*/

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

void H_to_bb_Virtual::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
}

void H_to_bb_Virtual::Calculate_ME2(const vector<Flavour>& flavs) {
  double value;
  double g_s = std::sqrt(4 * std::acos(-1) * alpha_qcd);  // strong gauge coupling; std::acos(-1) = pi
  double G_F = 1.16637886e-5; // Fermi constant in GeV^-2; value from PDG
  double vev = 1 / std::sqrt(G_F * std::sqrt(2)); // vacuum expectation value in GeV
  double mu = flavs[0].Mass(); // Renormalisation scale: Higgs mass in GeV
  double epsilon = 0.0001; // small parameter for dimensional regularization
  double m_b = flavs[2].Mass(); // b quark mass in GeV
  double constants = (-1) * std::pow(g_s * std::pow(mu, epsilon), 2) * m_b / vev;

}