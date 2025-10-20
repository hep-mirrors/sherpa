#include "EXTRA_XS/One2Three/H_to_bbg_Real.H"

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "EXTRA_XS/Main/ME_Tools.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include <assert.h>
#include "ATOOLS/Math/Vec4.H"
#include <cmath>
#include <numbers>

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace PHASIC;
using namespace std;

H_to_bbg_Real::H_to_bbg_Real(const vector<Flavour>& flavs, const Flavour& prop,
                     size_t non_prop, size_t gluon, size_t propj) :
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(4), m_anticur(4), m_nhel(4),
  m_prop(prop)
{
  DEBUG_FUNC(flavs<<" with prop "<<prop<<" in "<<gluon<<","<<propj);
  assert(non_prop>0 && gluon>0 && propj>0);
  if (flavs.size()!=4) THROW(fatal_error,"Internal error.");
  Vec4D k(1.0,0.0,1.0,0.0); // gauge

  for (size_t i(0);i<4;++i) { // iterate over the 4 flavours resp. particles
    // the Current_key object uniquely identify and retrieve a "current"; current = off-shell wave function/ building block to construct amplitude
    // it encapsulates information like flavour, charge, spin, polarization, etc.
    Current_Key ckey(flavs[i],MODEL::s_model,1);
    m_cur[i] = Current_Getter::GetObject("D"+ckey.Type(),ckey); // get the object that corresponds to the key "ckey"
    if (m_cur[i]==NULL) THROW(fatal_error, "current not found");
    m_cur[i]->SetDirection(i==0?1:-1); // +1 for i=0, -1 for i=1, 2, 3
    m_cur[i]->SetId(std::vector<int>(1,i));
    m_cur[i]->InitPols(std::vector<int>(1,m_spins[i])); // define allowed polarization states
    m_cur[i]->SetKey(i);
    m_cur[i]->SetGauge(k);
    m_nhel[i]=NHel(flavs[i]); // number of helicity states (based on spin properties)
  }
  // s-channel for prop (i,j)
  Current_Key ckey(prop,MODEL::s_model,2);
   m_scur = Current_Getter::GetObject("D"+ckey.Type(),ckey); // combine the two currents from i and j to one s-channel current
  METOOLS::Int_Vector isfs(2), ids(2), pols(2); // stores information about the outgoing particles: is fermion? Identifier, polarization
  isfs[0]=flavs[gluon].IsFermion();
  isfs[1]=flavs[propj].IsFermion();
  pols[0]=m_spins[ids[0]=gluon];
  pols[1]=m_spins[ids[1]=propj];
  m_scur->SetId(ids);
  m_scur->SetFId(isfs);
  m_scur->FindPermutations();
  // final current (1,2,3)
  ckey=Current_Key(flavs[0],MODEL::s_model,1);  // set up with incoming particle
  m_fcur = Current_Getter::GetObject("D"+ckey.Type(),ckey);
  METOOLS::Int_Vector isfs2(3), ids2(3), pols2(3); 
  isfs2[0]=flavs[1].IsFermion();
  isfs2[1]=flavs[2].IsFermion();
  isfs2[2]=flavs[3].IsFermion();
  pols2[0]=m_spins[ids2[0]=1]; // simultaneously sets up the ids
  pols2[1]=m_spins[ids2[1]=2];
  pols2[2]=m_spins[ids2[2]=3];
  m_fcur->SetId(ids2);
  m_fcur->SetFId(isfs2);
  m_fcur->FindPermutations(); // sets up allowed permutations of the outgoing particles
  // connect (2) & (3) into (2,3)
  m_v1=ConstructVertices(m_cur[gluon], m_cur[propj], m_scur); // vertex with b, gluon, bbar; bbar is actually an incoming b
  DEBUG_VAR(m_v1.size());
  // connect (1) & (2,3) into (1,2,3)
  m_v2=ConstructVertices(m_cur[non_prop],m_scur,m_fcur); //vertex with bbar, b, higgs
  DEBUG_VAR(m_v2.size());
  m_scur->Print();
  m_fcur->Print();
  m_scur->InitPols(pols);
  m_fcur->InitPols(pols2);
  m_fcur->HM().resize(m_n);
  for (size_t i(0);i<m_n;++i) m_fcur->HM()[i]=i;

  
  for (size_t i(0);i<4;++i) { // do the same for the anticurrent
    ckey=Current_Key(i==0?flavs[i]:flavs[i].Bar(),MODEL::s_model,1);
    m_anticur[i] = Current_Getter::GetObject("D"+ckey.Type(),ckey);
    if (m_anticur[i]==NULL) THROW(fatal_error, "current not found");
    m_anticur[i]->SetDirection(i==0?1:-1);
    m_anticur[i]->SetId(std::vector<int>(1,i));
    m_anticur[i]->InitPols(std::vector<int>(1,m_spins[i]));
    m_anticur[i]->SetKey(i);
    m_anticur[i]->SetGauge(k);
  }
  // s-channel for prop (2,3)
  ckey=Current_Key(prop.Bar(),MODEL::s_model,2);
  m_antiscur = Current_Getter::GetObject("D"+ckey.Type(),ckey);
  m_antiscur->SetId(ids);
  m_antiscur->SetFId(isfs);
  m_antiscur->FindPermutations();
  // final current (1,2,3)
  ckey=Current_Key(flavs[0].Bar(),MODEL::s_model,1);
  m_antifcur = Current_Getter::GetObject("D"+ckey.Type(),ckey);
  m_antifcur->SetId(ids2);
  m_antifcur->SetFId(isfs2);
  m_antifcur->FindPermutations();
  // connect (2) & (3) into (2,3)
  m_antiv1=ConstructVertices(m_anticur[gluon], m_anticur[propj], m_antiscur); // vertex with gluon, bbar, b; b is an incoming bbar
  DEBUG_VAR(m_antiv1.size());
  // connect (1) & (2,3) into (1,2,3)
  m_antiv2=ConstructVertices(m_anticur[non_prop],m_antiscur,m_antifcur); // vertex with b, bbar, higgs
  DEBUG_VAR(m_antiv2.size());
  m_antiscur->Print();
  m_antifcur->Print();
  m_antiscur->InitPols(pols);
  m_antifcur->InitPols(pols2);
  m_antifcur->HM().resize(m_n);
  for (size_t i(0);i<m_n;++i) m_antifcur->HM()[i]=i;

  p_ci=new Color_Integrator();
  Idx_Vector cids(4,0);
  METOOLS::Int_Vector acts(4,0), types(4,0); // (?) acts holds flags indicating whether each particle participates in strong interactions; 
  // (?) types stores the specific type of color charge for each particle
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


H_to_bbg_Real::~H_to_bbg_Real()
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


void H_to_bbg_Real::Born_setup(){}


void H_to_bbg_Real::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
  DEBUG_FUNC(momenta.size());
  // does not do anything yet because integrating this decay channel would result in infinities
  p_ci->GeneratePoint(); // create a new integration point for the color factors


  const std::vector<int> myI = { 0, 2, 1, 0 };
  const std::vector<int> myJ = { 0, 1, 0, 2 };

  //p_ci->SetI(myI);
  //p_ci->SetJ(myJ);

  if (anti) {
    for (size_t i(0);i<m_anticur.size();++i) {
      m_anticur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci->I()[i],p_ci->J()[i],0);
      m_anticur[i]->Print();
    }
    m_antiscur->Evaluate();
    m_antifcur->Evaluate();
  }
  else {
    for (size_t i(0);i<m_cur.size();++i) {
      m_cur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci->I()[i],p_ci->J()[i],0);
      m_cur[i]->Print();
    }
    m_scur->Evaluate();
    m_fcur->Evaluate();
  }

  vector<int> fill(m_n,1); // output amplitude vector
  for (size_t i(0);i<m_n;++i) (*this)[i]=Complex(0.0,0.0);
  if (anti) {
    m_antifcur->Contract<double>(*m_anticur.front(),fill,*this,0); // compute amplitude
  }
  else {
    m_fcur->Contract<double>(*m_cur.front(),fill,*this,0);
  }

  for (size_t i=0; i<size(); ++i) {
    (*this)[i] *= sqrt(p_ci->GlobalWeight()); // scale the final numerical result appropriately with the color factor
  }
  //std::cout << "GlobalWeight = " << p_ci->GlobalWeight() << std::endl;

/*
  for (size_t i = 0; i < p_ci->I().size(); ++i) {
      std::cout << "I[" << i << "] = " << p_ci->I()[i] << std::endl;
  }
  for (size_t i = 0; i < p_ci->J().size(); ++i) {
      std::cout << "J[" << i << "] = " << p_ci->J()[i] << std::endl;
  }
*/
Calculate_born_subtraction(momenta, anti);
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
  double m_i = std::sqrt(p_i * p_i);
  double m_j = std::sqrt(p_j * p_j);
  double m_k = std::sqrt(p_k * p_k);
  double mu_i = mu_n(p_i, p_j, p_k, m_i);
  double mu_j = mu_n(p_i, p_j, p_k, m_j);
  double mu_k = mu_n(p_i, p_j, p_k, m_k);
  double bracket = 2 * mu_k * mu_k + (1 - mu_i*mu_i - mu_j*mu_j - mu_k*mu_k)*(1 - y_ijk(p_i, p_j, p_k));
  double numerator = std::sqrt(bracket * bracket - 4 * mu_k * mu_k);
  double denominator = (1 - mu_i*mu_i - mu_j*mu_j - mu_k*mu_k)*(1 - y_ijk(p_i, p_j, p_k));
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
  if((p_i * p_i) < 0.0000000001){
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


static double V_ijk(ATOOLS::Vec4<double> p_i, ATOOLS::Vec4<double> p_j, ATOOLS::Vec4<double> p_k){
  double alpha_qcd = MODEL::s_model -> ScalarFunction("alpha_S", 125.25*125.25); // at Higgs scale
  
  #ifdef M_PI
    double pi = M_PI;
  #else
    const double pi = 3.14159265358979323846;
  #endif

  double C_F = 4.0 / 3.0;
  ATOOLS::Vec4<double> Q = p_i + p_j + p_k;
  double m2_Q = Q * Q; // Q mass squared

  double prefactor = 8 * pi * alpha_qcd * C_F;
  double first_summand = 2 / (1 - z_j_tilde(p_i, p_j, p_k) * (1 - y_ijk(p_i, p_j, p_k) ));
  double second_summand = nu_ijk_tilde(p_i, p_j, p_k) / nu_ijk(p_i, p_j, p_k) * (1 + z_j_tilde(p_i, p_j, p_k) + m2_Q / (p_i * p_j));

  return prefactor * (first_summand - second_summand);
}


void H_to_bbg_Real::Calculate_born_subtraction(const ATOOLS::Vec4D_Vector& momenta, bool anti){
  // implementation based on the formulas in the Catani Dittmaier Seymour Trocsanyi paper from 2002

  // first: some variables needed later (nomination also based on paper)
  ATOOLS::Vec4<double> p_g = momenta[1];
  ATOOLS::Vec4<double> p_b = momenta[2];
  ATOOLS::Vec4<double> p_bb = momenta[3];

  double V_gb_bb = V_ijk(p_g, p_b, p_bb);
  double V_gbb_b = V_ijk(p_g, p_bb, p_b);

  double born_ME2 = 9.1018603234124864; // value taken from the H -> bbb calculation in Comix1to2, value right out of Decay_Channel::ME2(...)
                                        // => not multiplied with any symmetry factors/ colour factors
  double m2_ij = p_b * p_b; // because m_i = 0 (gluon) and m_b = m_bb

  double D_gb_bb = V_gb_bb/ ((p_g + p_b)*(p_g + p_b) - m2_ij) * born_ME2;
  double D_gbb_b = V_gbb_b/ ((p_g + p_bb)*(p_g + p_bb) - m2_ij) * born_ME2;

  double subtraction_term = D_gb_bb + D_gbb_b;
}


size_t H_to_bbg_Real::NHel(const Flavour& fl)
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

