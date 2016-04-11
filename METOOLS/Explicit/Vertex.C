#include "METOOLS/Explicit/Vertex.H"

#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"

#include <algorithm>
#include <typeinfo>

using namespace METOOLS;
using namespace ATOOLS;

namespace METOOLS {

  template <class Type> std::string
  GetName(const Type &o,const int mode=0)
  {
    std::string id=Demangle(typeid(o).name());
    size_t pos=id.find("METOOLS::");
    if (pos<std::string::npos) id.erase(pos,9);
    pos=id.find("_Calculator");
    if (pos<std::string::npos) id.erase(pos,11);
    if (mode&1) {
      pos=id.find('<');
      if (pos<std::string::npos) {
	size_t epos=id.rfind('>');
	if (epos<std::string::npos) id.erase(pos,epos-pos+1);
      }
    }
    return id;
  }

}

size_t Vertex::s_vlmode(0);

std::map<std::string,Int_Vector> Vertex::s_h;

Vertex::Vertex(const Vertex_Key &key): 
  p_v(key.p_mv), p_c(NULL),
  p_info(key.p_dinfo), p_kin(NULL), p_h(NULL),
  m_sign(false), m_fperm(0), m_icplfac(1.0)
{
  if (key.p_mv==NULL) return;
  if (p_info)
    p_kin = new Dipole_Kinematics
      (p_info,key.m_j[0],key.m_j[1],key.p_k,key.p_c,key.p_kt);
  key.p_v=this;
  Vertex_Key ckey(key);
  for (ckey.m_n=0;ckey.m_n<key.p_mv->Lorentz.size();++ckey.m_n) {
    std::string ctag(ToString(ckey.p_mv->Color[ckey.m_n].PID()));
    if (key.p_dinfo) {
      if (abs(ckey.p_c->Flav().StrongCharge())==3) ctag="S-T";
      else if (key.p_c->Flav().StrongCharge()==8) ctag="S-F";
      else {
	for (size_t i(0);i<m_lc.size();++i) {
	  delete m_lc[i];
	  delete m_cc[i];
	}
	m_lc.clear();
	m_cc.clear();
	return;
      }
    }
    m_cc.push_back(CC_Getter::GetObject(ctag,ckey));
    if (m_cc.back()==NULL) {
      msg_Info()<<*ckey.p_mv<<std::endl;
      THROW(fatal_error,"Color calculator not implemented '"+
	    ctag+"'");
    }
    ckey.p_cc=m_cc.back();
    std::string lname=ckey.m_p;
    if(key.p_dinfo)
      lname+="X"+MODEL::s_model->MappedLorentzName(ckey.p_mv->Lorentz[ckey.m_n]);
    else lname+=ckey.p_mv->Lorentz[ckey.m_n];
    m_lc.push_back(LC_Getter::GetObject(lname,ckey));
    if (m_lc.back()==NULL) {
      msg_Out()<<*ckey.p_mv<<std::endl;
      THROW(fatal_error,"Lorentz calculator not implemented '"+lname+"'");
    }
  }
}

Vertex::~Vertex()
{
  for (size_t i(0);i<m_lc.size();++i) {
    delete m_lc[i];
    delete m_cc[i];
  }
  for (size_t i(0);i<m_j.size();++i)
    if (m_j[i]!=NULL) m_j[i]->DetachOut(this);
  if (p_kin) delete p_kin;
}

void Vertex::Evaluate()
{
  SetZero();
  if (p_kin && !p_kin->Trig()) return;
  for (Current_Vector::const_iterator jit(m_j.begin());
       jit!=m_j.end();++jit) if ((*jit)->Zero()) return;
  if (p_kin) {
    for (LC_Vector::const_iterator lit(m_lc.begin());
	 lit!=m_lc.end();++lit) (*lit)->Evaluate();
    if (!p_c->Zero()) {
      const CObject *c(p_kin->JK()->J().front().front());
      p_kin->JKT()->ConstructJ(p_kin->JKT()->P(),0,(*c)(0),(*c)(1),0);
    }
    return;
  }
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"():\n";
  msg_Indent();
#endif
  if (m_j.size()==2) {
    size_t hid(0), sh0(m_j[0]->J().size()), sh1(m_j[1]->J().size());
    for (size_t h0(0);h0<sh0;++h0) {
      const CObject_Vector *hjj0(&m_j[0]->J()[h0]);
      if (hjj0->empty()) { hid+=sh1; continue; }
      for (size_t h1(0);h1<sh1;++h1) {
	const CObject_Vector *hjj1(&m_j[1]->J()[h1]);
	if (hjj1->empty()) { ++hid; continue; }
	for (size_t c0(0);c0<hjj0->size();++c0) {
	  m_cjj[0]=(*hjj0)[c0];
	  for (size_t c1(0);c1<hjj1->size();++c1) {
	    m_cjj[1]=(*hjj1)[c1];
	    for (size_t k(0);k<m_cc.size();++k)
	      if (m_cc[k]->Evaluate(m_cjj)) {
		CObject *j(m_lc[k]->Evaluate(m_cjj));
		if (j==NULL) continue;
		j->Multiply(p_v->Coupling(k)*m_cc[k]->Coupling());
		j->SetH(H(hid));
		m_cc[k]->AddJ(j);
		SetZero(false);
	      }
	  }
	}
	++hid;
      }
    }
    return;
  }
  size_t hid(0);
  Int_Vector m_cjc(m_j.size()), m_hjc(m_j.size(),0);
  std::vector<const CObject_Vector*> m_hjj(m_j.size());
  for (size_t j(0);j<m_hjj.size();++j) m_hjj[j]=&m_j[j]->J().front();
  for (size_t hc(m_hjc.size()-1);m_hjc[0]<m_j[0]->J().size();) {
    if(m_hjc[hc]==m_j[hc]->J().size()){m_hjc[hc--]=0;++m_hjc[hc];continue;}
    m_hjj[hc]=&m_j[hc]->J()[m_hjc[hc]];if(hc<m_hjc.size()-1){++hc;continue;}
    for (Int_Vector::iterator i(m_cjc.begin());i!=m_cjc.end();++i) *i=0;
    bool zero(false);
    for (size_t i(0);i<m_cjj.size();++i)
      if (m_hjj[i]->empty()) {zero=true;break;}
      else m_cjj[i]=m_hjj[i]->front();
    if (zero) {++m_hjc[hc];++hid;continue;}
    for (size_t cc(m_cjc.size()-1);m_cjc[0]<m_hjj[0]->size();) {
      if (m_cjc[cc]==m_hjj[cc]->size()){m_cjc[cc--]=0;++m_cjc[cc];continue;}
      m_cjj[cc]=(*m_hjj[cc])[m_cjc[cc]];if(cc<m_cjc.size()-1){++cc;continue;}
      for (size_t k(0);k<m_cc.size();++k)
	if (m_cc[k]->Evaluate(m_cjj)) {
	  CObject *j(m_lc[k]->Evaluate(m_cjj));
	  if (j==NULL) continue;
	  j->Multiply(p_v->Coupling(k)*m_cc[k]->Coupling());
	  j->SetH(H(hid));
	  m_cc[k]->AddJ(j);
	  SetZero(false);
	}
      ++m_cjc[cc];
    }
    ++m_hjc[hc];
    ++hid;
  }
}

void Vertex::FindPermutation()
{
  m_fperm=0;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  Int_Vector id(p_c->Id()), fid(p_c->FId());
  Int_Vector pid(m_j[0]->Id()), pfid(m_j[0]->FId());
  for (size_t i(1);i<m_j.size();++i) {
    pid.insert(pid.end(),m_j[i]->Id().begin(),m_j[i]->Id().end());
    pfid.insert(pfid.end(),m_j[i]->FId().begin(),m_j[i]->FId().end());
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"  pid = "<<pid<<", pfid = "<<pfid<<"\n";
  msg_Debugging()<<"  id  = "<<id<<", fid  = "<<fid<<"\n";
#endif
  for (size_t i(0);i<id.size();++i) {
    for (size_t j(0);j<pid.size();++j) 
      if (pid[j]==id[i] && i!=j) {
	for (size_t k(j);k!=i;) {
	  size_t l(j>i?k-1:k+1);
 	  m_fperm+=pfid[k]==1&&pfid[l]==1;
#ifdef DEBUG__BG
	  if (pfid[k]==1 && pfid[l]==1)
	    msg_Debugging()<<"  swap "<<pid[l]<<" & "<<pid[k]<<"\n";
#endif
	  std::swap<int>(pid[k],pid[l]);
	  std::swap<int>(pfid[k],pfid[l]);
	  k=l;
	}
	break;
      }
  }
  m_sign=m_fperm%2==1;
#ifdef DEBUG__BG
  msg_Debugging()<<"} => "<<*this<<"\n";
#endif
}

void Vertex::InitPols()
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"() {\n";
#endif
  m_cjj.resize(m_j.size());
  int nmax(0);
  std::string id;
  for (size_t i(0);i<m_j.size();++i) {
    id+=ToString(m_j[i]->H().Spins());
    nmax=Max(nmax,m_j[i]->Id().back());
  }
  for (size_t i(0);i<=nmax;++i)
    for (size_t j(0);j<m_j.size();++j)
      if (std::find(m_j[j]->Id().begin(),
		    m_j[j]->Id().end(),i)!=
	  m_j[j]->Id().end()) {
	id+="_"+ToString(j);
	break;
      }
  std::map<std::string,Int_Vector>::iterator hit(s_h.find(id));
  if (hit!=s_h.end()) {
    p_h=&hit->second;
#ifdef DEBUG__BG
    msg_Debugging()<<"  "<<id<<" mapped to '"<<p_h<<"'\n}\n";
#endif
    return;
  }
  p_h=&s_h.insert(make_pair(id,Int_Vector())).first->second;
#ifdef DEBUG__BG
  msg_Debugging()<<"  "<<id<<" stored in '"<<p_h<<"'\n";
#endif
  Int_Vector m_hjc(m_j.size(),0);
  std::vector<Int_Vector> hjj(m_j.size());
  for (size_t i(0);i<hjj.size();++i) hjj[i]=m_j[i]->H()(0);
  for (size_t hc(m_hjc.size()-1);m_hjc[0]<m_j[0]->H().N();) {
    if(m_hjc[hc]==m_j[hc]->H().N()){m_hjc[hc--]=0;++m_hjc[hc];continue;}
    hjj[hc]=m_j[hc]->H()(m_hjc[hc]);if(hc<m_hjc.size()-1){++hc;continue;}
    Int_Vector ch(hjj.back()), id(m_j.back()->Id());
    for (size_t i(0);i<hjj.size()-1;++i) {
      for (size_t m(0);m<hjj[i].size();++m) {
	Int_Vector::iterator cit(ch.begin()), iit(id.begin());
	for (;iit<id.end();++iit,++cit) if (m_j[i]->Id()[m]<*iit) break;
	id.insert(iit,m_j[i]->Id()[m]);
	ch.insert(cit,hjj[i][m]);
      }
    }
#ifdef DEBUG__BG
    msg_Debugging()<<"  ["<<p_h->size()<<"]: j = "<<m_hjc
		   <<", h = "<<ch<<" -> id = "<<p_c->H()(ch)<<"\n";
#endif
    p_h->push_back(p_c->H()(ch));
    ++m_hjc[hc];
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}

bool Vertex::Map(const Vertex &v)
{
#ifdef DEBUG__BG
  msg_Debugging()<<"      "<<(m_cc.size()?GetName(*m_cc.front()):"")
		 <<"|"<<(m_lc.size()?GetName(*m_lc.front()):"")
		 <<" "<<VId()<<" "<<CVLabel()<<"\n"
		 <<"      "<<(v.m_cc.size()?GetName(*v.m_cc.front()):"")
		 <<"|"<<(v.m_lc.size()?GetName(*v.m_lc.front()):"")
		 <<" "<<v.VId()<<" "<<v.CVLabel()<<"\n";
#endif
  for (size_t i(0);i<m_cc.size();++i) {
    if (typeid(*m_cc[i])!=typeid(*v.m_cc[i])) return false;
    if (typeid(*m_lc[i])!=typeid(*v.m_lc[i])) return false;
    if (p_v->cpl[i].Value()!=v.p_v->cpl[i].Value()) return false;
  }
  return VId()==v.VId();
}

void Vertex::AddJ(const Current_Vector &j)
{
  for (size_t i(0);i<j.size();++i) AddJ(j[i]);
}

std::string Vertex::VId() const
{
  std::string estr("v");
  for (size_t i(0);i<m_j.size();++i)
    estr+="_"+ToString(m_j[i]->CId());
  return estr+"_"+ToString(p_c->CId());
}

std::string Vertex::CVLabel() const
{
  std::string label(m_lc[0]->Label());
  for (size_t i(1);i<m_lc.size();++i) label+=";"+m_lc[i]->Label();
  return label;
}

std::string Vertex::VLabel() const
{
  std::string label;
  if (s_vlmode&1) label+="\\scriptstyle\\blue F="+ToString(m_fperm);
  if (s_vlmode&2) {
    if (m_cc.empty() || m_lc.empty()) THROW(fatal_error,"Invalid call");
    std::string id(GetName(*m_cc.front())+"_"+GetName(*m_lc.front(),1));
    for (size_t pos;(pos=id.find("_"))!=std::string::npos && 
	   id[pos-1]!='\\';id.replace(pos,1,"\\_"));
    label+=std::string(label.length()>0?"\\\\":"")+
      "\\scriptstyle\\green T="+id+"("+m_j[0]->Flav().TexName();
    for (size_t i(1);i<m_j.size();++i) label+=","+m_j[i]->Flav().TexName();
    label+=")";
  }
  if (s_vlmode&4)
    label+=std::string(label.length()>0?"\\\\":"")+
      "\\scriptstyle\\red C="+CVLabel();
  for (size_t pos(label.find(','));
       pos!=std::string::npos;pos=label.find(',',pos+2))
    label.replace(pos,1,",,");
  return "decor.size=0ex,label=$\\begin{array}{c}"+label+"\\end{array}$";
}

void Vertex::CollectGraphs(Graph_Node *graph) const
{
  graph->push_back("    \\fmfv{"+VLabel()+"}{"+VId()+"}");
  graph->push_back("    %% "+VId());
  for (size_t i(0);i<m_j.size();++i) m_j[i]->CollectGraphs(graph);
}

const std::vector<int> &Vertex::Order() const
{
  return p_v->order;
}

int Vertex::Order(const size_t &id) const
{
  return p_v->order[id];
}

std::ostream &METOOLS::operator<<(std::ostream &str,const Vertex &v)
{
  for (size_t i(0);i<v.J().size();++i) {
    if (i) str<<"(+)";
    str<<'{'<<v.J(i)->Type()<<','<<v.J(i)->Flav()<<'}'<<v.J(i)->Id();
    if (v.J(i)->Sub()) str<<"S["<<v.J(i)->Sub()->Id()
			  <<v.J(i)->Sub()->Sub()->Id()<<"]";
  }
  if (v.JC()!=NULL) {
    str<<"-";
    if (v.Color().size() && v.Lorentz().size()) {
      str<<"'"<<GetName(*v.Color().front())
	 <<"*"<<GetName(*v.Lorentz().front());
      for (size_t i(1);i<v.Color().size();++i)
	str<<"+"<<GetName(*v.Color()[i])
	   <<"*"<<GetName(*v.Lorentz()[i]);
      str<<"'";
    }
    str<<v.Order()<<"->{"<<v.JC()->Type()
       <<','<<v.JC()->Flav()<<'}'<<v.JC()->Id();
  }
  if (v.Kin()) str<<" D["<<v.Kin()->JK()->Id()
		  <<","<<v.Kin()->Type()<<"]";
  return str<<" {"<<v.FPerm()<<","<<v.Sign()<<"}";
}
