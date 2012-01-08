#include "ATOOLS/Phys/Cluster_Amplitude.H"

#include <algorithm>
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Flow.H"

using namespace ATOOLS;

ClusterAmplitude_PVector::ClusterAmplitude_PVector()
{
#ifdef USING__Threading
  pthread_mutex_init(&m_mtx,NULL);
#endif
}

ClusterAmplitude_PVector::~ClusterAmplitude_PVector()
{
#ifdef USING__Threading
  pthread_mutex_destroy(&m_mtx);
#endif
  while (!empty()) {
    Cluster_Amplitude *ampl(back());
    pop_back();
    delete ampl;
  }
}

ClusterAmplitude_PVector Cluster_Amplitude::s_ampls;

Cluster_Amplitude::Cluster_Amplitude(Cluster_Amplitude *const prev):
  p_prev(prev), p_next(NULL), 
  m_oew(0), m_oqcd(0), m_nin(0), m_new(0), m_kin(0), m_nlo(0),
  m_mur2(0.0), m_muf2(0.0), m_mu2(0.0), m_kt2(0.0), m_z(0.0), m_phi(0.0),
  p_jf(NULL), p_procs(NULL), p_dinfo(NULL),
  p_rb(NULL), p_ms(NULL)
{
  if (p_prev!=NULL) p_prev->p_next=this;
}

Cluster_Amplitude::~Cluster_Amplitude()
{
  if (p_next) p_next->Delete();
  for (size_t i(0);i<m_legs.size();++i) m_legs[i]->Delete();
  if (p_prev) p_prev->p_next=NULL;
}

Cluster_Amplitude *Cluster_Amplitude::New
(Cluster_Amplitude *const prev)
{
  s_ampls.MtxLock();
  if (s_ampls.empty()) {
    s_ampls.MtxUnLock();
    return new Cluster_Amplitude(prev);
  }
  Cluster_Amplitude *ca(s_ampls.back());
  s_ampls.pop_back();
  s_ampls.MtxUnLock();
  ca->p_prev=prev;
  ca->p_next=NULL;
  ca->m_oew=ca->m_oqcd=0;
  ca->m_nin=ca->m_new=ca->m_kin=ca->m_nlo=0;
  ca->m_mur2=ca->m_muf2=ca->m_mu2=0.0;
  ca->m_kt2=ca->m_z=ca->m_phi=0.0;
  ca->p_jf=ca->p_procs=ca->p_dinfo=NULL;
  ca->p_rb=NULL;
  ca->p_ms=NULL;
  if (ca->p_prev!=NULL) ca->p_prev->p_next=ca;
  return ca;
}

void Cluster_Amplitude::Delete()
{
  if (p_next) p_next->Delete();
  for (size_t i(0);i<m_legs.size();++i) m_legs[i]->Delete();
  m_legs.clear();
  m_decs.clear();
  m_cmap.clear();
  if (p_prev) p_prev->p_next=NULL;
  p_prev=p_next=NULL;
  s_ampls.MtxLock();
  s_ampls.push_back(this);
  s_ampls.MtxUnLock();
}

void Cluster_Amplitude::CreateLeg
(const Vec4D &p,const Flavour &fl,
 const ColorID &col,const size_t &id)
{
  m_legs.push_back(Cluster_Leg::New(this,p,fl,col));
  if (id!=std::string::npos) m_legs.back()->SetId(id);
  else m_legs.back()->SetId(1<<(m_legs.size()-1));
}

Cluster_Amplitude *Cluster_Amplitude::Copy() const
{
  Cluster_Amplitude *copy(Cluster_Amplitude::New());
  copy->CopyFrom(this);
  return copy;
}

void Cluster_Amplitude::CopyFrom
(const Cluster_Amplitude *const master,const int mode)
{
  Cluster_Amplitude *sprev(p_prev), *snext(p_next);
  *this=*master;
  p_prev=sprev;
  p_next=snext;
  if (mode==1) m_legs.clear();
  else for (size_t i(0);i<m_legs.size();++i)
    m_legs[i] = Cluster_Leg::New(this,*master->m_legs[i]);
}

Cluster_Amplitude *Cluster_Amplitude::CopyNext() const
{
  const Cluster_Amplitude *root(this);
  Cluster_Amplitude *prev(NULL), *ref(NULL);
  while (root) {
    Cluster_Amplitude *copy(root->Copy());
    if (prev!=NULL) (prev->p_next=copy)->p_prev=prev;
    prev=copy;
    if (root==this) ref=prev;
    root=root->Next();
  }
  return ref;
}

Cluster_Amplitude *Cluster_Amplitude::CopyAll() const
{
  const Cluster_Amplitude *root(this);
  Cluster_Amplitude *prev(NULL), *ref(NULL);
  while (root->Prev()) root=root->Prev();
  while (root) {
    Cluster_Amplitude *copy(root->Copy());
    if (prev!=NULL) (prev->p_next=copy)->p_prev=prev;
    prev=copy;
    if (root==this) ref=prev;
    root=root->Next();
  }
  return ref;
}

void Cluster_Amplitude::CombineLegs
(Cluster_Leg *const i,Cluster_Leg *const j,
 const Flavour &fl,const ColorID &col)
{
  if (i->Amplitude()!=this || j->Amplitude()!=this) 
    THROW(fatal_error,"Leg not owned by current amplitude");
  for (ClusterLeg_Vector::iterator clit(m_legs.begin());
       clit!=m_legs.end();++clit) {
    if (*clit==i || *clit==j) {
      *clit = Cluster_Leg::New(this,i->Mom()+j->Mom(),fl,col);
      i->Delete();
      j->Delete();
      for (++clit;clit!=m_legs.end();++clit)
	if (*clit==i || *clit==j) {
	  clit=m_legs.erase(clit);
	  break;
	}
      break;
    }
  }
}

Cluster_Amplitude *Cluster_Amplitude::InitNext()
{
  if (p_next!=NULL) p_next->Delete();
  p_next = Cluster_Amplitude::New(this);
  return p_next;
}

Cluster_Amplitude *Cluster_Amplitude::InitPrev()
{
  if (p_prev!=NULL) return NULL;
  p_prev = new Cluster_Amplitude();
  p_prev->p_next=this;
  return p_prev;
}

void Cluster_Amplitude::SetNext(Cluster_Amplitude *const next) 
{
  if (p_next!=NULL) p_next->Delete();
  if (next->p_prev) next->p_prev->p_next=NULL;
  (p_next=next)->p_prev=this;
}

void Cluster_Amplitude::UnsetNext()
{
  if (p_next) p_next->p_prev=NULL;
  p_next=NULL;
}

void Cluster_Amplitude::DeletePrev()
{
  if (p_prev==NULL) return;
  p_prev->p_next=NULL;
  while (p_prev->Prev()) p_prev=p_prev->Prev();
  p_prev->Delete();
  p_prev=NULL;
}

void Cluster_Amplitude::DeleteNext()
{
  if (p_next!=NULL) p_next->Delete();
  p_next=NULL;
}

class Order_LegId {
public:
  int operator()(const Cluster_Leg *a,const Cluster_Leg *b) 
  { return a->Id()<b->Id(); }
};// end of class Order_LegId

void Cluster_Amplitude::IdSort()
{
  std::stable_sort(m_legs.begin(),m_legs.end(),Order_LegId());
}

size_t Cluster_Amplitude::NQCD() const
{
  size_t nqcd(0);
  for (size_t i(0);i<m_legs.size();++i)
    nqcd+=m_legs[i]->Flav().Strong();
  return nqcd;
}

size_t Cluster_Amplitude::NEW() const
{
  size_t nw(0);
  for (size_t i(0);i<m_legs.size();++i)
    nw+=!m_legs[i]->Flav().Strong();
  return nw;
}

Cluster_Leg *Cluster_Amplitude::IdLeg(const size_t &id) const
{
  for (size_t i(0);i<m_legs.size();i++)
    if (m_legs[i]->Id()==id) return m_legs[i];
  return NULL;
}

void Cluster_Amplitude::SetColours
(Cluster_Leg *const lij,Cluster_Leg *const li,Cluster_Leg *const lj)
{
  ColorID colij(lij->Col()), coli(0,0), colj(0,0);
  if (lij->Flav().StrongCharge()==3) {
    if (li->Flav().StrongCharge()==3) {
      if (lj->Flav().Strong()) {
	size_t nc(Flow::Counter());
	colj.m_j=coli.m_i=nc;
	colj.m_i=colij.m_i;
      }
      else {
	colj.m_j=colj.m_i=0;
	coli.m_i=colij.m_i;
      }
    }
    else {
      if (li->Flav().Strong()) {
	size_t nc(Flow::Counter());
	coli.m_j=colj.m_i=nc;
	coli.m_i=colij.m_i;
      }
      else {
	coli.m_j=coli.m_i=0;
	colj.m_i=colij.m_i;
      }
    }
  }
  else if (lij->Flav().StrongCharge()==-3) {
    if (li->Flav().StrongCharge()==-3) {
      if (lj->Flav().Strong()) {
	size_t nc(Flow::Counter());
	colj.m_i=coli.m_j=nc;
	colj.m_j=colij.m_j;
      }
      else {
	colj.m_j=colj.m_i=0;
	coli.m_j=colij.m_j;
      }
    }
    else {
      if (li->Flav().Strong()) {
	size_t nc(Flow::Counter());
	coli.m_i=colj.m_j=nc;
	coli.m_j=colij.m_j;
      }
      else {
	coli.m_j=coli.m_i=0;
	colj.m_j=colij.m_j;
      }
    }
  }
  else if (lij->Flav().Strong()) {
    if (li->Flav().StrongCharge()==8) {
      size_t nc(Flow::Counter());
      colj.m_i=coli.m_j=nc;
      colj.m_j=colij.m_j;
      coli.m_i=colij.m_i;
    }
    else if (abs(li->Flav().StrongCharge())==3) {
      coli.m_i=colij.m_i;
      colj.m_j=colij.m_j;
      if (li->Flav().StrongCharge()<0)
	std::swap<ColorID>(coli,colj);
    }
    else {
      coli.m_i=colij.m_i;
      coli.m_j=colij.m_j;
      if (lj->Flav().StrongCharge()==8)
	std::swap<ColorID>(coli,colj);
    }
  }
  else {
    if (li->Flav().Strong()) {
      size_t nc(Flow::Counter());
      coli.m_i=colj.m_j=nc;
      if (li->Flav().StrongCharge()<0)
	std::swap<ColorID>(coli,colj);
    }
  }
  li->SetCol(coli);
  lj->SetCol(colj);
}


size_t Cluster_Amplitude::IdIndex(const size_t &id) const
{
  for (size_t i(0);i<m_legs.size();i++)
    if (m_legs[i]->Id()==id) return i;
  return -1;
}

namespace ATOOLS {

  class Order_LegID {
  public:
    bool operator()(Cluster_Leg *const a,Cluster_Leg *const b) const
    { return a->Id()<b->Id(); }
  };// end of class Sort_LegID

}

void Cluster_Amplitude::OrderLegs()
{
  std::stable_sort(m_legs.begin(),m_legs.end(),Order_LegID());
}

namespace ATOOLS {

  std::ostream &operator<<
    (std::ostream &ostr,const Cluster_Amplitude &ampl)
  {
    ostr<<"("<<&ampl<<"): "<<ampl.NIn()
	<<" -> "<<ampl.Legs().size()-ampl.NIn()<<" {\n";
    ostr<<"  \\mu_r = "<<sqrt(ampl.MuR2())
	<<", \\mu_f = "<<sqrt(ampl.MuF2())
	<<", \\mu = "<<sqrt(ampl.Mu2())<<"\n";
    ostr<<"  k_T = "<<sqrt(ampl.KT2())<<", z = "<<ampl.Z()
	<<", phi = "<<ampl.Phi()<<", kin = "<<ampl.Kin()<<"\n";
    ostr<<"  oew = "<<ampl.OrderEW()<<", oqcd = "<<ampl.OrderQCD()
	<<", nlo = "<<ampl.NLO()<<", new = "<<ID(ampl.IdNew())<<"\n";
    if (ampl.Decays().size()) {
      std::string ds;
      for (DecayInfo_Vector::const_iterator cit(ampl.Decays().begin());
	   cit!=ampl.Decays().end();++cit) 
        ds+=ToString(**cit)+" ";
      ostr<<"  decs = { "<<ds<<"}\n";
    }
    if (ampl.RBMap()) {
      ostr<<"  rbs = { ";
      for (RB_Map::const_iterator rbit(ampl.RBMap()->begin());
	   rbit!=ampl.RBMap()->end();++rbit) 
	ostr<<"["<<rbit->first.m_type<<";"
	    <<rbit->first.m_flij<<"->"<<rbit->first.m_fli<<","
	    <<rbit->first.m_flj<<"]->("<<rbit->second->m_max<<") ";
      ostr<<"}\n";
    }
    if (ampl.ColorMap().size()) {
      std::string cs;
      for (CI_Map::const_iterator cit(ampl.ColorMap().begin());
	   cit!=ampl.ColorMap().end();++cit) 
	cs+=ToString(cit->first)+"->"+ToString(cit->second)+" ";
      ostr<<"  cols = { "<<cs<<"}\n";
    }
    for (size_t i(0);i<ampl.Legs().size();++i)
      ostr<<"  "<<*ampl.Legs()[i]<<"\n";
    return ostr<<"}";
  }

}
