#include "EXTRA_XS/One2Three/Comix1to3.H"

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include <assert.h>

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace PHASIC;
using namespace std;

Comix1to3::Comix1to3(const vector<Flavour>& flavs, const Flavour& prop,
                     size_t nonprop, size_t propi, size_t propj) :
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(4), m_anticur(4), m_nhel(4),
  m_prop(prop)
{
  DEBUG_FUNC(flavs<<" with prop "<<prop<<" in "<<propi<<","<<propj);
  assert(nonprop>0 && propi>0 && propj>0);
  if (flavs.size()!=4) THROW(fatal_error,"Internal error.");
  Vec4D k(1.0,0.0,1.0,0.0);

  for (size_t i(0);i<4;++i) {
    Current_Key ckey(i==0?flavs[i].Bar():flavs[i],MODEL::s_model,1);
    m_cur[i] = Current_Getter::GetObject("D"+ckey.Type(),ckey);
    if (m_cur[i]==NULL) THROW(fatal_error, "current not found");
    m_cur[i]->SetDirection(i==0?1:-1);
    m_cur[i]->SetId(std::vector<int>(1,i));
    m_cur[i]->InitPols(std::vector<int>(1,m_spins[i]));
    m_cur[i]->SetKey(i);
    m_cur[i]->SetGauge(k);
    m_nhel[i]=NHel(flavs[i]);
  }
  // s-channel for prop (i,j)
  Current_Key ckey(prop,MODEL::s_model,2);
  m_scur = Current_Getter::GetObject("D"+ckey.Type(),ckey);
  METOOLS::Int_Vector isfs(2), ids(2), pols(2);
  isfs[0]=flavs[propi].IsFermion();
  isfs[1]=flavs[propj].IsFermion();
  pols[0]=m_spins[ids[0]=propi];
  pols[1]=m_spins[ids[1]=propj];
  m_scur->SetId(ids);
  m_scur->SetFId(isfs);
  m_scur->FindPermutations();
  // final current (1,2,3)
  ckey=Current_Key(flavs[0],MODEL::s_model,1);
  m_fcur = Current_Getter::GetObject("D"+ckey.Type(),ckey);
  METOOLS::Int_Vector isfs2(3), ids2(3), pols2(3);
  isfs2[0]=flavs[1].IsFermion();
  isfs2[1]=flavs[2].IsFermion();
  isfs2[2]=flavs[3].IsFermion();
  pols2[0]=m_spins[ids2[0]=1];
  pols2[1]=m_spins[ids2[1]=2];
  pols2[2]=m_spins[ids2[2]=3];
  m_fcur->SetId(ids2);
  m_fcur->SetFId(isfs2);
  m_fcur->FindPermutations();
  // connect (2) & (3) into (2,3)
  m_v1=GetVertex(m_cur[propi], m_cur[propj], m_scur);
  DEBUG_VAR(*m_v1);
  // connect (1) & (2,3) into (1,2,3)
  m_v2=GetVertex(m_cur[nonprop],m_scur,m_fcur);
  DEBUG_VAR(*m_v2);
  m_scur->Print();
  m_fcur->Print();
  m_scur->InitPols(pols);
  m_fcur->InitPols(pols2);
  m_fcur->HM().resize(m_n);
  for (size_t i(0);i<m_n;++i) m_fcur->HM()[i]=i;

  for (size_t i(0);i<4;++i) {
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
  m_antiv1=GetVertex(m_anticur[propi], m_anticur[propj], m_antiscur);
  DEBUG_VAR(*m_antiv1);
  // connect (1) & (2,3) into (1,2,3)
  m_antiv2=GetVertex(m_anticur[nonprop],m_antiscur,m_antifcur);
  DEBUG_VAR(*m_antiv2);
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
    cids[i]=i;
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

Comix1to3::~Comix1to3()
{
  for (size_t i(0);i<4;++i) {
    delete m_cur[i];
    delete m_anticur[i];
  }
  delete m_scur;
  delete m_antiscur;
  delete m_fcur;
  delete m_antifcur;
}

Vertex* Comix1to3::GetVertex(Current* cur1, Current* cur2, Current* prop) {
  Vertex* v1(NULL);
  Current_Vector curs(2);
  curs[0]=cur1;
  curs[1]=cur2;
  Vertex_Key vkey(curs,prop,MODEL::s_model);
  MODEL::VMIterator_Pair keyrange(MODEL::s_model->GetVertex(vkey.ID()));
  DEBUG_VAR(vkey.ID());
  if (keyrange.first!=keyrange.second) {
    vkey.p_mv=keyrange.first->second;//fixme?
    vkey.m_p=std::string(1,'D');
    v1=new Vertex(vkey);
  }
  else {
    std::swap<Current*>(curs[0],curs[1]);
    vkey=Vertex_Key(curs,prop,MODEL::s_model);
    DEBUG_VAR(vkey.ID());
    keyrange=MODEL::s_model->GetVertex(vkey.ID());
    if (keyrange.first!=keyrange.second) {
      vkey.p_mv=keyrange.first->second;//fixme?
      vkey.m_p=std::string(1,'D');
      v1=new Vertex(vkey);
    }
    else THROW(fatal_error, "vertex not found: "+vkey.ID());
  }
  v1->AddJ(vkey.m_j);
  v1->SetJC(prop);
  return v1;
}
      

void Comix1to3::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
  DEBUG_FUNC(momenta.size());
  p_ci->GeneratePoint();
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


  vector<int> fill(m_n,1);
  for (size_t i(0);i<m_n;++i) (*this)[i]=Complex(0.0,0.0);
  if (anti) {
    m_antifcur->Contract<double>(*m_anticur.front(),fill,*this,0);
  }
  else {
    m_fcur->Contract<double>(*m_cur.front(),fill,*this,0);
  }

  for (size_t i=0; i<size(); ++i) {
    (*this)[i] *= sqrt(p_ci->GlobalWeight());
  }
}

size_t Comix1to3::NHel(const Flavour& fl)
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
