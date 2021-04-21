#include "PHASIC++/Selectors/Selector.H"

namespace ATOOLS {

  template <typename NType> class Variable_Base;

  class Order_Base;

}

namespace PHASIC {

  class Variable_Selector : public Selector_Base {

    ATOOLS::Variable_Base<double>    *p_variable;
    std::vector<ATOOLS::Order_Base*>  m_orders;

    std::vector<std::pair<double,double> > m_bounds;

    std::vector<std::vector<ATOOLS::Vec4D_Vector> > m_moms;
    std::vector<ATOOLS::Flavour_Vector> m_cfl;
    std::vector<std::vector<size_t> > m_nfl;
    std::vector<size_t> m_ffl;

    std::string m_omode;

    int m_imode;

    bool Trigger(const ATOOLS::Selector_List &sl,size_t &l,size_t &u,
		 ATOOLS::Vec4D_Vector &moms,const size_t &f,
		 const size_t &n,const size_t &m,const size_t &id);

  public:

    // constructor
    Variable_Selector(Process_Base* const proc,
                      const int& imode,
                      const std::string& name);

    // destructor
    ~Variable_Selector();

    // member functions
    void BuildCuts(Cut_Data *cuts);

    void SetRange(int id,ATOOLS::Flavour_Vector pfl,
		  ATOOLS::Flavour_Vector fl,
		  std::vector<std::pair<double,double> > &bounds);

    bool Trigger(ATOOLS::Selector_List &sl);
    bool Trigger(ATOOLS::Selector_List &sl,const int id);

  };// end of class Variable_Selector

}// end of namespace PHASIC

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Math/Variable.H"
#include "ATOOLS/Phys/Ordering.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <algorithm>

#define DEBUG__Variable_Selector

using namespace PHASIC;
using namespace ATOOLS;

Variable_Selector::Variable_Selector(Process_Base* const proc,
                                     const int& imode,
                                     const std::string& name) :
  Selector_Base("Variable("+name+")",proc)
{
  m_imode=imode;
  size_t pos(name.find('|'));
  while (pos<name.length()-1 && name[pos+1]=='|') pos=name.find('|',pos+2);
  std::string vname(name.substr(0,pos));
  p_variable = ATOOLS::Variable_Getter::GetObject(vname,vname);
  if (p_variable==NULL) THROW
    (fatal_error,"Variable '"+vname+"' does not exist. Run 'Sherpa"+
       std::string(" SHOW_VARIABLE_SYNTAX=1' to list variables."));
  m_omode=name.substr(pos!=std::string::npos?pos+1:pos);
  if (m_omode!="") {
    if (m_omode[0]=='[') {
      if (m_omode[m_omode.length()-1]!=']') 
	THROW(fatal_error,"Invalid ordering mode '"+m_omode+"'");
      Data_Reader reader(",",":","!","=");
      std::string mode(m_omode.substr(1));
      mode.erase(mode.length()-1,1);
      if (mode.length()>0) {
	reader.SetString(mode);
	std::vector<std::string> omodes;
	if (!reader.VectorFromString(omodes,""))
	  THROW(critical_error,"Invalid ordering mode '"+m_omode+"'");
	for (size_t i(0);i<omodes.size();++i) {
	  m_orders.push_back(Order_Getter::GetObject(omodes[i],""));
	  if (m_orders.back()==NULL) 
	    THROW(fatal_error,"Invalid ordering mode '"+omodes[i]+"'");
	}
      }
    }
    else if (m_omode[0]=='{') {
      if (m_omode[m_omode.length()-1]!='}') 
	THROW(fatal_error,"Invalid ordering mode '"+m_omode+"'");
      Data_Reader reader(",",":","!","=");
      std::string ffl(m_omode.substr(1));
      ffl.erase(ffl.length()-1,1);
      if (ffl.length()>0) {
	reader.SetString(ffl);
	if (!reader.VectorFromString(m_ffl,""))
	  THROW(critical_error,
		"Invalid Syntax in Selector.dat: '"+m_omode+"'");
      }
    }
  }
}

Variable_Selector::~Variable_Selector() 
{
  while (m_orders.size()) {
    delete m_orders.back();
    m_orders.pop_back();
  }
  delete p_variable;
}

void Variable_Selector::BuildCuts(Cut_Data *cuts)
{
}

void Variable_Selector::SetRange
(int id,std::vector<Flavour> pfl,std::vector<Flavour> fl,
 std::vector<std::pair<double,double> > &bounds)
{
  m_cfl.push_back(Flavour_Vector());
  m_nfl.push_back(std::vector<size_t>());
  if (id>=m_cfl.size()) THROW(fatal_error,"Invalid call");
  for (size_t i(0);i<fl.size();++i) {
    bool found(false);
    for (size_t j(0);j<m_cfl[id].size();++j)
      if (m_cfl[id][j]==fl[i]) {
	++m_nfl[id][j];
	found=true;
	break;
      }
    if (!found) {
      m_cfl[id].push_back(fl[i]);
      m_nfl[id].push_back(1);
    }
  }
  m_moms.push_back(std::vector<Vec4D_Vector>(m_cfl[id].size()));
  m_bounds=bounds;
  m_name="Variable_Selector_"+ToString(m_imode)+"_";
  for (size_t j(0);j<m_cfl[id].size();++j) {
    m_name+="_"+m_cfl[id][j].IDName()+"-"+ToString(m_nfl[id][j]);
    for (int i(m_imode?0:m_nin);i<pfl.size();++i)
      if (m_cfl[id][j].Includes(pfl[i]))
	m_moms[id][j].push_back(Vec4D());
  }
  msg_Debugging()<<METHOD<<"(): order = "<<m_omode
		 <<", imode = "<<m_imode<<" {\n";
  for (size_t j(0);j<m_bounds.size();++j) {
    msg_Debugging()<<"  "<<p_variable->Name()<<"_{"<<j<<"}";
    if (m_ffl.size()>j) msg_Debugging()<<"["<<m_ffl[j]<<"]";
    if (m_orders.size()>j) msg_Debugging()<<"["<<m_orders[j]<<"]";
    msg_Debugging()<<" -> "<<m_bounds[j].first
		   <<" .. "<<m_bounds[j].second<<"\n";
  }
  for (size_t j(0);j<m_cfl[id].size();++j) {
    msg_Debugging()<<"  "<<j<<": "<<m_cfl[id][j].IDName()
		   <<" ("<<m_nfl[id][j]<<") -> {";
    if (m_moms[id][j].size()>0) msg_Debugging()<<m_moms[id][j].front();
    for (size_t k(1);k<m_moms[id][j].size();++k)
      msg_Debugging()<<","<<m_moms[id][j][k];
    msg_Debugging()<<"}\n";
  }
  msg_Debugging()<<"}\n";
  m_sel_log->ChangeName(m_name);
}

bool Variable_Selector::Trigger
(const Selector_List &sl,size_t &l,size_t &u,std::vector<Vec4D> &moms,
 const size_t &f,const size_t &n,const size_t &m,const size_t &id) 
{
  msg_Indent();
  if (f==m_cfl[id].size()) {
    if (m_ffl.empty()) u=l;
    else if (u>=m_ffl.size() || l!=m_ffl[u]) {
      ++l;
      return true;
    }
    if (u>=m_bounds.size()) return true;
    double v((*p_variable)(&moms.front(),moms.size()));
#ifdef DEBUG__Variable_Selector
    msg_Debugging()<<u<<"th ("<<l<<") "<<p_variable->Name()
		   <<"="<<v<<" vs. {"<<m_bounds[u].first
		   <<","<<m_bounds[u].second<<"}\n";
#endif
    bool res(v<m_bounds[u].first || v>m_bounds[u].second);
    ++l; ++u;
    return !m_sel_log->Hit(res);
  }
  if (n==m_nfl[id][f]) return Trigger(sl,l,u,moms,f+1,0,0,id);
  moms.push_back(Vec4D());
  for (size_t k(m);k<m_moms[id][f].size();++k) {
#ifdef DEBUG__Variable_Selector
    msg_Debugging()<<"f = "<<f<<", n = "<<n<<", m = "<<m
		   <<", k = "<<k<<" -> "<<m_cfl[id][f].IDName()
		   <<" ("<<m_moms[id][f][k]<<") {\n";
#endif
    moms.back()=m_moms[id][f][k];
    if (!Trigger(sl,l,u,moms,f,n+1,k+1,id)) return false;
#ifdef DEBUG__Variable_Selector
    msg_Debugging()<<"}\n";
#endif
  }
  moms.pop_back();
  return true;
}

bool Variable_Selector::Trigger(Selector_List &sl,const int id)
{
#ifdef DEBUG__Variable_Selector
  msg_Debugging()<<METHOD<<"(id="<<0<<"): {\n";
#endif
  for (size_t j(0);j<m_cfl[id].size();++j) {
    size_t i(0);
    for (size_t k(m_imode?0:m_nin);k<sl.size();++k)
      if (m_cfl[id][j].Includes(sl[k].Flavour()))
	m_moms[id][j][i++]=sl[k].Momentum();
    while (i<m_moms[id][j].size()) m_moms[id][j][i++]=Vec4D();
  }
  size_t l(0), u(0);
  std::vector<Vec4D> moms;
  bool hit(Trigger(sl,l,u,moms,0,0,0,0));
#ifdef DEBUG__Variable_Selector
  msg_Debugging()<<"}\n";
#endif
  return hit;
}

bool Variable_Selector::Trigger(Selector_List &sl)
{
  return Trigger(sl,p_sub?(p_sub->IsReal()?0:p_sub->m_idx+1):0);
}

DECLARE_GETTER(Variable_Selector,"VariableSelector",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Variable_Selector>::
operator()(const Selector_Key &key) const
{
#ifdef DEBUG__Variable_Selector
  msg_Debugging()<<"Getter<Variable_Selector>::operator(): {\n";
#endif

  auto s = key.m_settings["VariableSelector"];
  s.DeclareVectorSettingsWithEmptyDefault({ "Flavs" });
  s.DeclareMatrixSettingsWithEmptyDefault({ "Ranges" });

  const auto flavs = s["Flavs"].GetVector<int>();
  if (flavs.empty())
    THROW(critical_error,"Missing \"Flav\" specification in variable selector");
  Flavour_Vector cflavs(flavs.size());
  for (size_t j(0);j<flavs.size();++j) {
    cflavs[j]=Flavour((kf_code)abs(flavs[j]));
    if (flavs[j]<0) cflavs[j]=cflavs[j].Bar();
  }

  const auto bounds = s["Ranges"].GetMatrix<double>();
  if (bounds.empty())
    THROW(critical_error,"Missing \"Ranges\" specification in variable selector");
  std::vector<std::pair<double,double> > cbounds;
  for (const auto single_bounds : bounds) {
    if (single_bounds.size() != 2)
      THROW(critical_error,"Ranges need to have two entries.");
    cbounds.push_back(std::make_pair(single_bounds[0], single_bounds[1]));
  }

  auto name = s["Variable"].SetDefault("").Get<std::string>();
  const auto orderings = s["Ordering"].SetDefault("").Get<std::string>();
  name += "|" + orderings;
  const auto imode = s["Mode"].SetDefault(0).Get<int>();
  Variable_Selector *vs(
      new Variable_Selector(key.p_proc, imode, name));
  vs->SetRange(0,key.p_proc->Flavours(),cflavs,cbounds);
  NLO_subevtlist *subs(key.p_proc->GetSubevtList());
  if (subs) {
    for (size_t i(0);i<subs->size()-1;++i) {
      Flavour_Vector fls((*subs)[i]->p_fl,
			 &(*subs)[i]->p_fl[(*subs)[i]->m_n]);
      vs->SetRange((*subs)[i]->m_idx+1,fls,cflavs,cbounds);
    }
  }
#ifdef DEBUG__Variable_Selector
  msg_Debugging()<<"}\n";
#endif
  return vs;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Variable_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"variable selector";
}
