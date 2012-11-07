#include "EXTRA_XS/One2Two/Comix1to2.H"

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Comix1to2::Comix1to2(const vector<Flavour>& flavs) :
  Spin_Amplitudes(flavs,Complex(0.0,0.0)), m_cur(3), m_anticur(3), m_nhel(3)
{
  if (flavs.size()!=3) THROW(fatal_error,"Internal error.");

  Vec4D k(1.0,0.0,1.0,0.0);

  for (size_t i(0);i<3;++i) {
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
  // final current (1,2)
  Current_Key ckey(flavs[0],MODEL::s_model,1);
  m_fcur = Current_Getter::GetObject("D"+ckey.Type(),ckey);
  Int_Vector isfs(2), ids(2), pols(2);
  isfs[0]=flavs[1].IsFermion();
  isfs[1]=flavs[2].IsFermion();
  pols[0]=m_spins[ids[0]=1];
  pols[1]=m_spins[ids[1]=2];
  m_fcur->SetId(ids);
  m_fcur->SetFId(isfs);
  m_fcur->FindPermutations();
  // connect (1) & (2) into (1,2)
  m_v1=GetVertex(m_cur[1], m_cur[2], m_fcur);
  m_fcur->InitPols(pols);
  m_fcur->Print();
  m_fcur->HM().resize(m_n);
  for (size_t i(0);i<m_n;++i) m_fcur->HM()[i]=i;

  for (size_t i(0);i<3;++i) {
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
  m_antiv1=GetVertex(m_anticur[1], m_anticur[2], m_antifcur);
  m_antifcur->InitPols(pols);
  m_antifcur->Print();
  m_antifcur->HM().resize(m_n);
  for (size_t i(0);i<m_n;++i) m_antifcur->HM()[i]=i;
}


Comix1to2::~Comix1to2()
{
  for (size_t i(0);i<3;++i) {
    delete m_cur[i];
    delete m_anticur[i];
  }
  delete m_fcur;
  delete m_antifcur;
}

Vertex* Comix1to2::GetVertex(Current* cur1, Current* cur2, Current* prop) {
  Vertex* v1(NULL);
  Vertex_Key vkey(cur1,cur2,NULL,prop,MODEL::s_model);
  MODEL::VMIterator_Pair keyrange(MODEL::s_model->GetVertex(vkey.ID()));
  if (keyrange.first!=keyrange.second && keyrange.first->second->on) {
    vkey.p_mv=keyrange.first->second;//fixme!!
    vkey.m_p=std::string(1,'D');
    v1=new Vertex(vkey);
  }
  else {
    vkey=Vertex_Key(cur2,cur1,NULL,prop,MODEL::s_model);
    keyrange=MODEL::s_model->GetVertex(vkey.ID());
    if (keyrange.first!=keyrange.second && keyrange.first->second->on) {
      vkey.p_mv=MODEL::s_model->GetVertex(vkey.ID()).first->second;//fixme!!
      vkey.m_p=std::string(1,'D');
      v1=new Vertex(vkey);
    }
    else THROW(fatal_error, "vertex not found: "+vkey.ID());
  }
  v1->SetJA(vkey.p_a);
  v1->SetJB(vkey.p_b);
  v1->SetJC(prop);
  return v1;
}
     

void Comix1to2::Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti) {
  if (anti) {
    for (size_t i(0);i<m_anticur.size();++i) {
      m_anticur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,0,0);
      m_anticur[i]->Print();
    }
    m_antifcur->Evaluate();
  }
  else {
    for (size_t i(0);i<m_cur.size();++i) {
      m_cur[i]->ConstructJ(i==0?-momenta[i]:momenta[i],0,0,0);
      m_cur[i]->Print();
    }
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

long unsigned int Comix1to2::IHM(const vector<int>& l)
{
  long unsigned int ihm=0;
  for (size_t i=0; i<l.size(); ++i) {
    int li=l[i];
    switch (m_nhel[i]) {
    case 1:
      break;
    case 2:
        if (li==0) ihm+=1<<i;
      break;
    case 3:
        if (li==0 || li==2) ihm+=1<<i;
      break;
    default:
      THROW(not_implemented, "Not yet capable of spin > 1.");
    }
  }
  return ihm;
}

long unsigned int Comix1to2::IHP(const vector<int>& l)
{
  long unsigned int ihp=0;
  for (size_t i=0; i<l.size(); ++i) {
    int li=l[i];
    switch (m_nhel[i]) {
    case 1:
      break;
    case 2:
        if (li==1) ihp+=1<<i;
      break;
    case 3:
        if (li==1 || li==2) ihp+=1<<i;
      break;
    default:
      THROW(not_implemented, "Not yet capable of spin > 1.");
    }
  }
  return ihp;
}
