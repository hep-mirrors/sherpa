#include "EXTRA_XS/One2Two/Comix1to2.H"

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "EXTRA_XS/Main/ME_Tools.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "EXTRA_XS/One2Two/H_to_bb_Virtual.H"
#include <memory>
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace PHASIC;
using namespace std;

Comix1to2::Comix1to2(const vector<Flavour>& flavs) :
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(3), m_anticur(3), m_nhel(3), m_flavs(flavs)
{
  if (flavs.size()!=3) THROW(fatal_error,"Internal error.");

  IsNLODecay(); // sets the isNLO flag true if this decay is an NLO decay. False otherwise.

  Vec4D k(1.0,0.0,1.0,0.0); // gauge

  for (size_t i(0);i<3;++i) { // iterate over the 3 flavours
    // the Current_key object uniquely identify and retrieve a "current"; current = off-shell wave function/ building block to construct amplitude
    // it encapsulates information like flavour, charge, spin, polarization, etc.
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
  // final current (1,2): represents the combined state of the two outgoing particles
  Current_Key ckey(flavs[0],MODEL::s_model,1); // set up with incoming particle
  m_fcur = Current_Getter::GetObject("D"+ckey.Type(),ckey);
  METOOLS::Int_Vector isfs(2), ids(2), pols(2); // stores information about the outgoing particles: is fermion? Identifier, polarization
  isfs[0]=flavs[1].IsFermion();
  isfs[1]=flavs[2].IsFermion();
  pols[0]=m_spins[ids[0]=1]; // simultaneously sets up the ids
  pols[1]=m_spins[ids[1]=2];
  m_fcur->SetId(ids);
  m_fcur->SetFId(isfs);
  m_fcur->FindPermutations(); // sets up allowed permutations of the outgoing particles
  // connect (1) & (2) into (1,2)
  m_v1=ConstructVertices(m_cur[1], m_cur[2], m_fcur);
  DEBUG_VAR(m_v1.size());
  m_fcur->InitPols(pols);
  m_fcur->Print();
  m_fcur->HM().resize(m_n);
  for (size_t i(0);i<m_n;++i) m_fcur->HM()[i]=i;

  for (size_t i(0);i<3;++i) { // do the same for the anticurrent
    ckey=Current_Key(i==0?flavs[i]:flavs[i].Bar(),MODEL::s_model,1);
    m_anticur[i] = Current_Getter::GetObject("D"+ckey.Type(),ckey);
    if (m_anticur[i]==NULL) THROW(fatal_error, "current not found");
    m_anticur[i]->SetDirection(i==0?1:-1);
    m_anticur[i]->SetId(std::vector<int>(1,i));
    m_anticur[i]->InitPols(std::vector<int>(1,m_spins[i]));
    m_anticur[i]->SetKey(i);
    m_anticur[i]->SetGauge(k);
  }
  // final current (1,2)
  ckey=Current_Key(flavs[0].Bar(),MODEL::s_model,1);
  m_antifcur = Current_Getter::GetObject("D"+ckey.Type(),ckey);
  m_antifcur->SetId(ids);
  m_antifcur->SetFId(isfs);
  m_antifcur->FindPermutations();
  // connect (1) & (2) into (1,2)
  m_antiv1=ConstructVertices(m_anticur[1], m_anticur[2], m_antifcur);
  DEBUG_VAR(m_antiv1.size());
  m_antifcur->InitPols(pols);
  m_antifcur->Print();
  m_antifcur->HM().resize(m_n);
  for (size_t i(0);i<m_n;++i) m_antifcur->HM()[i]=i;

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


Comix1to2::~Comix1to2()
{
  for (size_t i(0);i<3;++i) {
    delete m_cur[i];
    delete m_anticur[i];
  }
  delete m_fcur;
  delete m_antifcur;
  if (p_ci) delete p_ci;
}


bool Comix1to2::IsNLODecay() {
  // check if the decay channel is decaying at NLO
  // so far only h0 -> b bbar at NLO is included
  if (m_flavs[0].IDName() == "h0" && m_flavs[1].IDName() == "b" && m_flavs[2].IDName() == "bb"){
    isNLO = true;
  } else{
    isNLO = false;
  }
  return isNLO;
}


void Comix1to2::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
  p_ci->GeneratePoint();
  if (anti) {
    for (size_t i(0);i<m_anticur.size();++i) {
      m_anticur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci->I()[i],p_ci->J()[i],0);
      m_anticur[i]->Print();
    }
    m_antifcur->Evaluate();
  }
  else {
    for (size_t i(0);i<m_cur.size();++i) {
      m_cur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,p_ci->I()[i],p_ci->J()[i],0);
      m_cur[i]->Print();
    }
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
}


double Comix1to2::get_NLO_part(){
  return 0;
}


METOOLS::Amplitude2_Tensor Comix1to2::AddNLOTensor(METOOLS::Amplitude2_Tensor old_tensor){
  return old_tensor;
}


std::string Comix1to2::getType(){
  return "LO";
}


size_t Comix1to2::NHel(const Flavour& fl)
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
    THROW(not_implemented, "Not yet capable of spin > 1.");
    return 0;
  }
}
