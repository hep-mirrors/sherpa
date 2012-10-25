#ifndef PHASIC_Selectors_Variable_Selector_H
#define PHASIC_Selectors_Variable_Selector_H

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

    std::vector<std::vector<int> >    m_sels;
    std::vector<ATOOLS::Vec4D_Vector> m_moms;
    ATOOLS::Flavour_Vector            m_cfl;
    std::vector<size_t>               m_nfl, m_ffl;

    std::string m_omode;

    int m_imode;

    bool Trigger(const ATOOLS::Vec4D_Vector &p,size_t &l,size_t &u,
		 ATOOLS::Vec4D_Vector &moms,const size_t &f,
		 const size_t &n,const size_t &m);

  public:

    // constructor
    Variable_Selector(const int &nin,const int &nout,
		      const int &imode,ATOOLS::Flavour *const fl,
		      const std::string &name);

    // destructor
    ~Variable_Selector();

    // member functions
    void BuildCuts(Cut_Data *cuts);

    void SetRange(ATOOLS::Flavour_Vector fl,
		  std::vector<std::pair<double,double> > &bounds);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);
    bool NoJetTrigger(const ATOOLS::Vec4D_Vector &p);
    bool JetTrigger(const ATOOLS::Vec4D_Vector &,
		    ATOOLS::NLO_subevtlist *const);

  };// end of class Variable_Selector

}// end of namespace PHASIC

#endif

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

Variable_Selector::Variable_Selector
(const int &nin,const int &nout,const int &imode,Flavour *const fl,
 const std::string &name): Selector_Base("Variable("+name+")")
{
  m_nin=nin;
  m_nout=nout;
  m_n=m_nin+m_nout;
  m_imode=imode;
  m_fl = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i]=fl[i];
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
  delete [] m_fl;
}

void Variable_Selector::BuildCuts(Cut_Data *cuts)
{
}

void Variable_Selector::SetRange
(std::vector<Flavour> fl,std::vector<std::pair<double,double> > &bounds)
{
  for (size_t i(0);i<fl.size();++i) {
    bool found(false);
    for (size_t j(0);j<m_cfl.size();++j)
      if (m_cfl[j]==fl[i]) {
	++m_nfl[j];
	found=true;
	break;
      }
    if (!found) {
      m_cfl.push_back(fl[i]);
      m_nfl.push_back(1);
    }
  }
  m_sels.resize(m_cfl.size());
  m_moms.resize(m_cfl.size());
  m_bounds=bounds;
  m_name="Variable_Selector_"+ToString(m_imode)+"_";
  for (size_t j(0);j<m_cfl.size();++j) {
    m_name+="_"+m_cfl[j].IDName()+"-"+ToString(m_nfl[j]);
    for (int i(m_imode?0:m_nin);i<m_n;++i)
      if (m_cfl[j].Includes(m_fl[i]))
	m_sels[j].push_back(i);
    m_moms[j].resize(m_sels[j].size());
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
  for (size_t j(0);j<m_cfl.size();++j) {
    msg_Debugging()<<"  "<<j<<": "<<m_cfl[j].IDName()
		   <<" ("<<m_nfl[j]<<") -> {";
    if (m_sels[j].size()>0) msg_Debugging()<<m_sels[j].front();
    for (size_t k(1);k<m_sels[j].size();++k)
      msg_Debugging()<<","<<m_sels[j][k];
    msg_Debugging()<<"}\n";
  }
  msg_Debugging()<<"}\n";
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool Variable_Selector::Trigger
(const Vec4D_Vector &p,size_t &l,size_t &u,std::vector<Vec4D> &moms,
 const size_t &f,const size_t &n,const size_t &m) 
{
  msg_Indent();
  if (f==m_cfl.size()) {
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
  if (n==m_nfl[f]) return Trigger(p,l,u,moms,f+1,0,0);
  moms.push_back(Vec4D());
  for (size_t k(m);k<m_sels[f].size();++k) {
#ifdef DEBUG__Variable_Selector
    msg_Debugging()<<"f = "<<f<<", n = "<<n<<", m = "<<m
		   <<", k = "<<k<<" -> "<<m_cfl[f].IDName()
		   <<" ("<<m_sels[f][k]<<") {\n";
#endif
    moms.back()=m_moms[f][k];
    if (!Trigger(p,l,u,moms,f,n+1,k+1)) return false;
#ifdef DEBUG__Variable_Selector
    msg_Debugging()<<"}\n";
#endif
  }
  moms.pop_back();
  return true;
}

bool Variable_Selector::Trigger(const Vec4D_Vector &p) 
{
#ifdef DEBUG__Variable_Selector
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  for (size_t j(0);j<m_cfl.size();++j) {
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
    if (m_orders.size()>j)
      std::sort(m_moms[j].begin(),m_moms[j].end(),*m_orders[j]);
  }
  size_t l(0), u(0);
  std::vector<Vec4D> moms;
  bool hit(Trigger(p,l,u,moms,0,0,0));
#ifdef DEBUG__Variable_Selector
  msg_Debugging()<<"}\n";
#endif
  return hit;
}

namespace PHASIC{

DECLARE_ND_GETTER(Variable_Selector_Getter,"\"",Selector_Base,Selector_Key,true);

Selector_Base *Variable_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<2) THROW(critical_error,"Invalid syntax");
  Data_Reader reader(",",":","!","=");
  reader.SetString(key[0][0]);
  std::vector<int> flavs;
  if (!reader.VectorFromString(flavs,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+key[0][0]+"'");
  Flavour_Vector cflavs(flavs.size());
  for (size_t j(0);j<flavs.size();++j) {
    cflavs[j]=Flavour((kf_code)abs(flavs[j]));
    if (flavs[j]<0) cflavs[j]=cflavs[j].Bar();
  }
  reader.SetString(key[0][1]);
  std::vector<std::vector<double> > crits;
  if (!reader.MatrixFromString(crits,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+key[0][1]+"'");
  std::vector<std::pair<double,double> > bounds;
  for (size_t j(0);j<crits.size();++j) {
    if (crits[j].size()<2) 
      THROW(critical_error,
	    "Invalid Syntax in Selector.dat: '"+key[0][1]+"'");
    bounds.push_back(std::pair<double,double>
		     (crits[j][0],crits[j][1]));
  }
  std::string tag(key.m_key.substr(1));
  tag.erase(tag.length()-1,1);
  tag+="|"+(key[0].size()>2?key[0][2]:"");
  int imode(0);
  if (key.front().size()>3) imode=ToType<int>(key[0][3]);
  Variable_Selector *vs(new Variable_Selector
			(key.p_proc->NIn(),key.p_proc->NOut(),imode,
			 (Flavour*)&key.p_proc->Process()->
			 Flavours().front(),tag));
  vs->SetRange(cflavs,bounds);
  vs->SetProcess(key.p_proc);
  return vs;
}

bool Variable_Selector::NoJetTrigger(const Vec4D_Vector &p) 
{
  return true;
}

bool Variable_Selector::JetTrigger
(const Vec4D_Vector &p,NLO_subevtlist *const)
{
  return Trigger(p);
}

void Variable_Selector_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"variable selector"; 
}

}
