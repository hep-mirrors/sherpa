#include "COMIX/Amplitude/Vertex_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE COMIX::Vertex_Key
#define OBJECT_TYPE COMIX::Vertex_Base
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

template class Getter_Function
<void,const COMIX::Model*,std::less<std::string> >;

using namespace COMIX;
using namespace ATOOLS;

size_t Vertex_Base::s_vlmode(0);
size_t Vertex_Base::s_cimin(1);
size_t Vertex_Base::s_cimax(3);

Vertex_Base::Vertex_Base(const Vertex_Key &key,
			 const size_t &oew,const size_t &oqcd): 
  p_a(NULL), p_b(NULL), p_c(NULL), 
  m_sign(false), m_act(true), 
  m_fperm(0), m_oew(oew), m_oqcd(oqcd),
  m_tag(key.Type()), m_cplfac(1.0) {}

Vertex_Base::~Vertex_Base()
{
  if (p_a!=NULL) p_a->DetachOut(this);
  if (p_b!=NULL) p_b->DetachOut(this);
}

void Vertex_Base::FindPermutation()
{
  m_fperm=0;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  Int_Vector id(p_c->Id()), fid(p_c->FId());
  Int_Vector pid(p_a->Id()), pfid(p_a->FId());
  if (1) {//p_b->Id().front()>id.front()) {
    pid.insert(pid.end(),p_b->Id().begin(),p_b->Id().end());
    pfid.insert(pfid.end(),p_b->FId().begin(),p_b->FId().end());
  }
  else {
    pid.insert(pid.begin(),p_b->Id().begin(),p_b->Id().end());
    pfid.insert(pfid.begin(),p_b->FId().begin(),p_b->FId().end());
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

bool Vertex_Base::Map(const Vertex_Base &v)
{
#ifdef DEBUG__BG
  msg_Debugging()<<"      "<<typeid(*this).name()
		 <<" "<<VId()<<" "<<CVLabel()<<"\n"
		 <<"      "<<typeid(v).name()
		 <<" "<<v.VId()<<" "<<v.CVLabel()<<"\n";
#endif
  if (typeid(*this)!=typeid(v)) return false;
  if (VId()!=v.VId()) return false;
  return CVLabel()==v.CVLabel();
}

std::string Vertex_Base::VId() const
{
  return "v_"+ToString(p_a->CId())+"_"
    +ToString(p_b->CId())+"_"+ToString(p_c->CId());
}

std::string Vertex_Base::CVLabel() const
{
  return "";
}

std::string Vertex_Base::VLabel() const
{
  std::string label;
  if (s_vlmode&1) label+="\\scriptstyle\\blue F="+ToString(m_fperm);
  if (s_vlmode&2) {
    std::string id(typeid(*this).name());
    size_t pos(id.find("COMIX"));
    if (pos<id.length()-7) id=id.substr(pos+6,id.length()-7-pos);
    while ((pos=id.find("_"))!=std::string::npos && 
	   id[pos-1]!='\\') id.replace(pos,1,"\\_");
    label+=std::string(label.length()>0?"\\\\":"")+
      "\\scriptstyle\\green T="+id+"("+p_a->Flav().TexName()+",,"+
      p_b->Flav().TexName()+",,"+p_c->Flav().TexName()+")";
  }
  if (s_vlmode&4)
    label+=std::string(label.length()>0?"\\\\":"")+
      "\\scriptstyle\\red C="+CVLabel();
  return "decor.size=0ex,label=$\\begin{array}{c}"+label+"\\end{array}$";
}

void Vertex_Base::CollectGraphs(Graph_Node *graph) const
{
  graph->push_back("    \\fmfv{"+VLabel()+"}{"+VId()+"}");
  graph->push_back("    %% "+VId());
  p_a->CollectGraphs(graph);
  p_b->CollectGraphs(graph);
}

std::string Vertex_Key::Type() const
{
  return std::string(1,p_a->Type())+p_b->Type()+p_c->Type();
}

std::string Vertex_Key::ID(const int dir) const
{
  Flavour a(p_a->Flav()), b(p_b->Flav()), c(p_c->Flav());
  switch(dir) {
  case 2:
    std::swap<Flavour>(a,b);
  case 1: 
    std::swap<Flavour>(b,c);
    b=b.Bar();
    c=c.Bar();
    break;
  case 0: break;
  default:
    THROW(fatal_error,"Internal error");
  }
  return '{'+a.IDName()+"}{"+b.IDName()+"}{"+c.IDName()+'}';
}

bool Vertex_Key::operator<(const Vertex_Key &k) const
{
  if (p_a<k.p_a) return true;
  if (p_a>k.p_a) return false;
  if (p_b<k.p_b) return true;
  if (p_b>k.p_b) return false;
  return false;
}

std::ostream &COMIX::operator<<(std::ostream &str,const Vertex_Base &v)
{
  if (v.JA()!=NULL) 
    str<<'{'<<v.JA()->Type()<<','<<v.JA()->Flav()<<'}'<<v.JA()->Id();
  if (v.JB()!=NULL) 
    str<<"(+){"<<v.JB()->Type()<<','<<v.JB()->Flav()<<'}'<<v.JB()->Id();
  if (v.JC()!=NULL) str<<"-'"<<typeid(v).name()<<"'("<<v.OrderEW()<<","
		       <<v.OrderQCD()<<")->{"<<v.JC()->Type()
		       <<','<<v.JC()->Flav()<<'}'<<v.JC()->Id();
  return str<<" {"<<v.FPerm()<<","<<v.Sign()<<"}";
}
