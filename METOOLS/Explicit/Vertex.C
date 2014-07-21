#include "METOOLS/Explicit/Vertex.H"

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
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

Vertex::Vertex(const Vertex_Key &key): 
  p_a(NULL), p_b(NULL), p_c(NULL), p_e(NULL),
  p_info(key.p_dinfo), p_kin(NULL),
  m_sign(false), m_act(true), 
  m_fperm(0), m_oew(0), m_oqcd(0),
  m_icplfac(1.0)
{
  if (key.p_mv==NULL) return;
  if (p_info)
    p_kin = new Dipole_Kinematics
      (p_info,key.p_a,key.p_b,key.p_k,key.p_c,key.p_kt);
  key.p_v=this;
  m_oew=key.p_mv->oew;
  m_oqcd=key.p_mv->oqcd;
  m_cpl.resize(key.p_mv->cpl.size());
  for (size_t i(0);i<key.p_mv->cpl.size();++i) 
    m_cpl[i]=key.p_mv->cpl[i].Value();
  Vertex_Key ckey(key);
  for (ckey.m_n=0;ckey.m_n<key.p_mv->Lorentz.size();++ckey.m_n) {
    std::string ctag(ToString(ckey.p_mv->Color[ckey.m_n].PID()));
    if (key.p_dinfo) {
      if (abs(ckey.p_c->Flav().StrongCharge())==3) ctag="S-T";
      else if (key.p_c->Flav().StrongCharge()==8) ctag="S-F";
      else {
	m_act=false;
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
    m_lc.push_back(LC_Getter::GetObject
		   (ckey.m_p+ckey.p_mv->Lorentz[ckey.m_n]->Type(),ckey));
    if (m_lc.back()==NULL) {
      msg_Out()<<*ckey.p_mv<<std::endl;
      THROW(fatal_error,"Lorentz calculator not implemented '"+
	    ckey.m_p+ckey.p_mv->Lorentz[ckey.m_n]->Type()+"'");
    }
  }
}

Vertex::~Vertex()
{
  for (size_t i(0);i<m_lc.size();++i) {
    delete m_lc[i];
    delete m_cc[i];
  }
  if (p_a!=NULL) p_a->DetachOut(this);
  if (p_b!=NULL) p_b->DetachOut(this);
  if (p_e!=NULL) p_e->DetachOut(this);
  if (p_kin) delete p_kin;
}

Current *Vertex::J(const size_t &i) const
{
  switch (i) {
  case 0: return p_a;
  case 1: return p_b;
  case 2: return p_e;
  default: THROW(fatal_error,"Invalid index "+ToString(i));
  }
  return NULL;
}

void Vertex::Evaluate()
{
  if (p_kin && !p_kin->Trig()) return;
  for (LC_Vector::const_iterator lit(m_lc.begin());
       lit!=m_lc.end();++lit) (*lit)->Evaluate();
  if (p_kin && !p_c->Zero()) {
    const CObject *c(p_kin->JK()->J().front().front());
    p_kin->JKT()->ConstructJ(p_kin->JKT()->P(),0,(*c)(0),(*c)(1),0);
  }
}

void Vertex::FindPermutation()
{
  m_fperm=0;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): swap = "
		 <<m_lc.front()->Swap()<<" {\n";
#endif
  Current *ja(p_a), *jb(p_b);
  if (m_lc.front()->Swap()) std::swap<Current*>(ja,jb);
  Int_Vector id(p_c->Id()), fid(p_c->FId());
  Int_Vector pid(ja->Id()), pfid(ja->FId());
  pid.insert(pid.end(),jb->Id().begin(),jb->Id().end());
  pfid.insert(pfid.end(),jb->FId().begin(),jb->FId().end());
  if (p_e!=NULL) {
    pid.insert(pid.end(),p_e->Id().begin(),p_e->Id().end());
    pfid.insert(pfid.end(),p_e->FId().begin(),p_e->FId().end());
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
  m_h.clear();
  if (p_e) {
    for (size_t i(0);i<p_a->H().N();++i) {
      const Int_Vector &ca(p_a->H()(i));
      for (size_t j(0);j<p_b->H().N();++j) {
	const Int_Vector cb(p_b->H()(j));
	for (size_t k(0);k<p_e->H().N();++k) {
	  Int_Vector ch(p_e->H()(k)), id(p_e->Id());
	  for (size_t m(0);m<p_b->Id().size();++m) {
	    Int_Vector::iterator cit(ch.begin()), iit(id.begin());
	    for (;iit<id.end();++iit,++cit) if (p_b->Id()[m]<*iit) break;
	    id.insert(iit,p_b->Id()[m]);
	    ch.insert(cit,cb[m]);
	  }
	  for (size_t m(0);m<p_a->Id().size();++m) {
	    Int_Vector::iterator cit(ch.begin()), iit(id.begin());
	    for (;iit<id.end();++iit,++cit) if (p_a->Id()[m]<*iit) break;
	    id.insert(iit,p_a->Id()[m]);
	    ch.insert(cit,ca[m]);
	  }
#ifdef DEBUG__BG
	  msg_Debugging()<<"  "<<ch<<" -> "<<p_c->H()(ch)<<"\n";
#endif
	  m_h.push_back(p_c->H()(ch));
	}
      }
    }
  }
  else {
    for (size_t i(0);i<p_a->H().N();++i) {
      const Int_Vector &ca(p_a->H()(i));
      for (size_t j(0);j<p_b->H().N();++j) {
	Int_Vector ch(p_b->H()(j)), id(p_b->Id());
	for (size_t m(0);m<p_a->Id().size();++m) {
	  Int_Vector::iterator cit(ch.begin()), iit(id.begin());
	  for (;iit<id.end();++iit,++cit) if (p_a->Id()[m]<*iit) break;
	  id.insert(iit,p_a->Id()[m]);
	  ch.insert(cit,ca[m]);
	}
#ifdef DEBUG__BG
	msg_Debugging()<<"  "<<ch<<" -> "<<p_c->H()(ch)<<"\n";
#endif
	m_h.push_back(p_c->H()(ch));
      }
    }
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
  for (size_t i(0);i<m_cc.size();++i)
    if (typeid(*m_cc[i])!=typeid(*v.m_cc[i])) return false;
  for (size_t i(0);i<m_lc.size();++i)
    if (typeid(*m_lc[i])!=typeid(*v.m_lc[i])) return false;
  if (VId()!=v.VId()) return false;
  return CVLabel()==v.CVLabel();
}

std::string Vertex::VId() const
{
  std::string estr;
  if (p_e!=NULL) estr="_"+ToString(p_e->CId());
  return "v_"+ToString(p_a->CId())+"_"
    +ToString(p_b->CId())+estr+"_"+ToString(p_c->CId());
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
    std::string estr;
    if (p_e!=NULL) estr=","+p_e->Flav().TexName();
    label+=std::string(label.length()>0?"\\\\":"")+
      "\\scriptstyle\\green T="+id+"("+p_a->Flav().TexName()+","+
      p_b->Flav().TexName()+estr+","+p_c->Flav().TexName()+")";
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
  p_a->CollectGraphs(graph);
  p_b->CollectGraphs(graph);
  if (p_e!=NULL) p_e->CollectGraphs(graph);
}

std::ostream &METOOLS::operator<<(std::ostream &str,const Vertex &v)
{
  if (v.JA()!=NULL) {
    str<<'{'<<v.JA()->Type()<<','<<v.JA()->Flav()<<'}'<<v.JA()->Id();
    if (v.JA()->Sub()) str<<"S["<<v.JA()->Sub()->Id()
			  <<v.JA()->Sub()->Sub()->Id()<<"]";
  }
  if (v.JB()!=NULL) {
    str<<"(+){"<<v.JB()->Type()<<','<<v.JB()->Flav()<<'}'<<v.JB()->Id();
    if (v.JB()->Sub()) str<<"S["<<v.JB()->Sub()->Id()
			  <<v.JB()->Sub()->Sub()->Id()<<"]";
  }
  if (v.JE()!=NULL) {
    str<<"(+){"<<v.JE()->Type()<<','<<v.JE()->Flav()<<'}'<<v.JE()->Id();
    if (v.JE()->Sub()) str<<"S["<<v.JE()->Sub()->Id()
			  <<v.JE()->Sub()->Sub()->Id()<<"]";
  }
  if (v.JC()!=NULL) {
    str<<"-";
    if (v.Color().size() && v.Lorentz().size()) str<<"'"<<
      GetName(*v.Color().front())<<"|"<<GetName(*v.Lorentz().front())<<"'";
    str<<"("<<v.OrderEW()<<","<<v.OrderQCD()<<")->{"
       <<v.JC()->Type()<<','<<v.JC()->Flav()<<'}'<<v.JC()->Id();
  }
  if (v.Kin()) str<<" D["<<v.Kin()->JK()->Id()
		  <<","<<v.Kin()->Type()<<"]";
  return str<<" {"<<v.FPerm()<<","<<v.Sign()<<"}";
}
