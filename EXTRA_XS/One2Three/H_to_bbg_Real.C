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
#include "EXTRA_XS/One2Three/CF_Decl.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"


using namespace EXTRAXS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace PHASIC;
using namespace std;

H_to_bbg_Real::H_to_bbg_Real(const vector<Flavour>& flavs, const Flavour& prop,
                     size_t non_prop, size_t gluon, size_t propj) :
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(4), m_anticur(4), m_nhel(4),
  m_prop(prop), m_non_prop(non_prop), m_gluon(gluon), m_propj(propj)
{
  DEBUG_FUNC(flavs<<" with prop "<<prop<<" in "<<m_gluon<<","<<m_propj);
  assert(m_non_prop>0 && m_gluon>0 && m_propj>0);
  if (flavs.size()!=4) THROW(fatal_error,"Internal error.");
  Vec4D k(1.0,0.0,1.0,0.0); // gauge

  scale = flavs[0].Mass(); // scale = mass of incoming particle
  CalculateAlphaQCD(scale);
  
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
  isfs[0]=flavs[m_gluon].IsFermion();
  isfs[1]=flavs[m_propj].IsFermion();
  pols[0]=m_spins[ids[0]=m_gluon];
  pols[1]=m_spins[ids[1]=m_propj];
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
  m_v1=ConstructVertices(m_cur[m_gluon], m_cur[m_propj], m_scur); // vertex with b, gluon, bbar; bbar is actually an incoming b
  DEBUG_VAR(m_v1.size());
  // connect (1) & (2,3) into (1,2,3)
  m_v2=ConstructVertices(m_cur[m_non_prop],m_scur,m_fcur); //vertex with bbar, b, higgs
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
  m_antiv1=ConstructVertices(m_anticur[m_gluon], m_anticur[m_propj], m_antiscur); // vertex with gluon, bbar, b; b is an incoming bbar
  DEBUG_VAR(m_antiv1.size());
  // connect (1) & (2,3) into (1,2,3)
  m_antiv2=ConstructVertices(m_anticur[m_non_prop],m_antiscur,m_antifcur); // vertex with b, bbar, higgs
  DEBUG_VAR(m_antiv2.size());
  m_antiscur->Print();
  m_antifcur->Print();
  m_antiscur->InitPols(pols);
  m_antifcur->InitPols(pols2);
  m_antifcur->HM().resize(m_n);
  for (size_t i(0);i<m_n;++i) m_antifcur->HM()[i]=i;

  p_ci=new Color_Integrator();
  Idx_Vector cids(4,0);
  METOOLS::Int_Vector acts(4,0), types(4,0);
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


bool H_to_bbg_Real::IsNLODecay(){
  return true;
}


void H_to_bbg_Real::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
  DEBUG_FUNC(momenta.size());
  p_ci->GeneratePoint(); // create a new integration point for the color factors

  const std::vector<int> myI = { 0, 2, 1, 0 };
  const std::vector<int> myJ = { 0, 1, 0, 2 };

  p_ci->SetI(myI);
  p_ci->SetJ(myJ);

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

  // falsely the m_z scale is used for the calculation of alpha_S instead of the Higgs scale. To correct this, multiply the calculated amplitude tensor (*this)
  // with \sqrt{alpha_S(M_H) / alpha_S(scale used here)}.

  double alpha_qcd_used = MODEL::s_model -> ScalarConstant("alpha_S"); // currently used alpha_S value
  double global_factor = alpha_qcd / alpha_qcd_used;
  for (size_t i=0; i<size(); ++i) {
   (*this)[i] *= std::sqrt(global_factor);
  }
}


double H_to_bbg_Real::get_NLO_ME2(){
  return 0;   // even though this is technically a NLO class, the real contribution is calculated analog to LO in the Decay_Channel.
              // Therefore, adding it's ME2 contribution via the get_NLO_part method is not necessary.
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


void H_to_bbg_Real::CalculateAlphaQCD(double scale) {
  alpha_qcd = (MODEL::s_model) -> ScalarFunction("alpha_S", scale*scale); // at Higgs scale
}


METOOLS::Amplitude2_Tensor H_to_bbg_Real::AddNLOTensor(METOOLS::Amplitude2_Tensor old_tensor){
  return old_tensor;
}


std::string H_to_bbg_Real::getType(){
  return "R";
}
