#include "ATOOLS/Phys/Cluster_Amplitude.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;

Cluster_Amplitude::Cluster_Amplitude(Cluster_Amplitude *const prev):
  p_prev(prev), p_next(NULL), 
  m_oew(0), m_oqcd(0), m_swap(0), m_nin(0), m_new(0),
  m_mur2(0.0), m_muf2(0.0), m_kt2qcd(0.0), 
  m_x1(1.0), m_x2(1.0), m_rbmax(1.0),
  p_jf(NULL)
{
  if (p_prev!=NULL) p_prev->p_next=this;
}

Cluster_Amplitude::~Cluster_Amplitude()
{
  if (p_next) delete p_next;
  for (size_t i(0);i<m_legs.size();++i) delete m_legs[i];
  if (p_prev) p_prev->p_next=NULL;
}

void Cluster_Amplitude::CreateLeg
(const Vec4D &p,const Flavour &fl,
 const ColorID &col,const size_t &id)
{
  m_legs.push_back(new Cluster_Leg(this,p,fl,col));
  if (id!=std::string::npos) m_legs.back()->SetId(id);
  else m_legs.back()->SetId(1<<(m_legs.size()-1));
}

Cluster_Amplitude *Cluster_Amplitude::Copy() const
{
  Cluster_Amplitude *copy(new Cluster_Amplitude());
  copy->CopyFrom(this);
  return copy;
}

void Cluster_Amplitude::CopyFrom
(const Cluster_Amplitude *const master)
{
  Cluster_Amplitude *sprev(p_prev), *snext(p_next);
  *this=*master;
  p_prev=sprev;
  p_next=snext;
  for (size_t i(0);i<m_legs.size();++i)
    m_legs[i] = new Cluster_Leg(this,*master->m_legs[i]);
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
      *clit = new Cluster_Leg(this,i->Mom()+j->Mom(),fl,col);
      delete i;
      delete j;
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
  if (p_next!=NULL) delete p_next;
  p_next = new Cluster_Amplitude(this);
  return p_next;
}

void Cluster_Amplitude::SetNext(Cluster_Amplitude *const next) 
{
  if (p_next!=NULL) delete p_next;
  if (next->p_prev) next->p_prev->p_next=NULL;
  (p_next=next)->p_prev=this;
}

void Cluster_Amplitude::DeletePrev()
{
  if (p_prev==NULL) return;
  p_prev->p_next=NULL;
  while (p_prev->Prev()) p_prev=p_prev->Prev();
  delete p_prev;
  p_prev=NULL;
}

void Cluster_Amplitude::DeleteNext()
{
  if (p_next!=NULL) delete p_next;
  p_next=NULL;
}

void Cluster_Amplitude::Print() const
{
  msg_Out()<<"("<<this<<")["<<m_swap<<"]: "<<m_nin
	   <<" -> "<<m_legs.size()-m_nin<<" {\n";
  msg_Out()<<"  \\mu_r = "<<sqrt(m_mur2)
	   <<", \\mu_f = "<<sqrt(m_muf2)<<"\n";
  msg_Out()<<"  x_1 = "<<m_x1<<", x_2 = "<<m_x2<<"\n";
  msg_Out()<<"  k_{T,QCD} = "<<sqrt(m_kt2qcd)
	   <<", (R/B)_{max} = "<<m_rbmax<<"\n";
  msg_Out()<<"  oew = "<<m_oew
	   <<", oqcd = "<<m_oqcd<<"\n";
  msg_Out()<<"  swap = "<<m_swap
	   <<", new = "<<ID(m_new)<<"\n";
  if (m_cmap.size()) {
    std::string cs;
    for (CI_Map::const_iterator cit(m_cmap.begin());
	 cit!=m_cmap.end();++cit) 
      cs+=ToString(cit->first)+"->"+ToString(cit->second)+" ";
    msg_Out()<<"  cols = { "<<cs<<"}\n";
  }
  for (size_t i(0);i<m_legs.size();++i)
    msg_Out()<<"  "<<*m_legs[i]<<"\n";
  msg_Out()<<"}\n";
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

void Cluster_Amplitude::SwapInOrder() 
{ 
  std::swap<Cluster_Leg*>(m_legs[0],m_legs[1]);
  std::swap<double>(m_x1,m_x2);
  m_swap=1-m_swap;
}

Cluster_Leg *Cluster_Amplitude::IdLeg(const size_t &id) const
{
  for (size_t i(0);i<m_legs.size();i++)
    if (m_legs[i]->Id()==id) return m_legs[i];
  return NULL;
}

size_t Cluster_Amplitude::IdIndex(const size_t &id) const
{
  for (size_t i(0);i<m_legs.size();i++)
    if (m_legs[i]->Id()==id) return i;
  return -1;
}

namespace ATOOLS {

  std::ostream &operator<<
    (std::ostream &ostr,const Cluster_Amplitude &ampl)
  {
    ostr<<"("<<&ampl<<")["<<ampl.InSwaped()<<"]: "<<ampl.NIn()
	<<" -> "<<ampl.Legs().size()-ampl.NIn()<<" {\n";
    ostr<<"  \\mu_r = "<<sqrt(ampl.MuR2())
	<<", \\mu_f = "<<sqrt(ampl.MuF2())<<"\n";
    ostr<<"  x_1 = "<<ampl.X1()<<", x_2 = "<<ampl.X2()<<"\n";
    ostr<<"  k_{T,QCD} = "<<sqrt(ampl.KT2QCD())
	<<", (R/B)_{max} = "<<ampl.RBMax()<<"\n";
    ostr<<"  oew = "<<ampl.OrderEW()
	<<", oqcd = "<<ampl.OrderQCD()<<"\n";
    msg_Out()<<"  swap = "<<ampl.InSwaped()
	     <<", new = "<<ID(ampl.New())<<"\n";
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
