#include "Variable_Selector.H"

#include "Variable.H"
#include "Ordering.H"
#include "MyStrStream.H"
#include "Exception.H"

#define DEBUG__Variable_Selector

using namespace ATOOLS;

Variable_Selector::Variable_Selector
(const int &nin,const int &nout,Flavour *const fl,
 const std::string &name): p_order(NULL)
{
  m_name="Variable("+name+")";
  m_nin=nin;
  m_nout=nout;
  m_n=m_nin+m_nout;
  m_fl = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i]=fl[i];
  std::string vname(name.substr(0,name.find('|')));
  p_variable = Variable_Getter::GetObject(vname,vname);
  if (p_variable==NULL) THROW
    (fatal_error,"Variable '"+vname+"' does not exist. Run 'Sherpa"+
       std::string(" SHOW_VARIABLE_SYNTAX=1' to list variables."));
  m_omode=name.substr(name.find('|')+1);
  if (m_omode!="") {
    p_order = Order_Getter::GetObject(m_omode,"");
    if (p_order==NULL) 
      THROW(fatal_error,"Invalid ordering mode '"+m_omode+"'");
  }
}

Variable_Selector::~Variable_Selector() 
{
  if (p_order!=NULL) delete p_order;
  delete p_variable;
  delete [] m_fl;
}

void Variable_Selector::BuildCuts(Cut_Data *cuts)
{
}

void Variable_Selector::SetRange
(std::vector<Flavour> fl,std::vector<std::pair<double,double> > &bounds)
{
  if (fl.size()<2) THROW(fatal_error,"Wrong number of flavours");
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
  m_name="Variable_Selector";
  for (size_t j(0);j<m_cfl.size();++j) {
    m_name+="_"+m_cfl[j].IDName()+"-"+ToString(m_nfl[j]);
    for (int i(m_nin);i<m_n;++i)
      if (m_cfl[j].Includes(m_fl[i])) m_sels[j].push_back(i);
    m_moms[j].resize(m_sels[j].size());
  }
  msg_Debugging()<<METHOD<<"(): order = "<<m_omode<<" {\n";
  for (size_t j(0);j<m_bounds.size();++j)
    msg_Debugging()<<"  "<<p_variable->Name()<<"_{"<<j
		   <<"} -> "<<m_bounds[j].first
		   <<" .. "<<m_bounds[j].second<<"\n";
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
(const Vec4D *p,size_t &l,std::vector<Vec4D> &moms,
 const size_t &f,const size_t &n,const size_t &m) 
{
  msg_Indent();
  if (f==m_cfl.size()) {
    if (l>=m_bounds.size()) return true;
    double v((*p_variable)(&moms.front(),moms.size()));
#ifdef DEBUG__Variable_Selector
    msg_Debugging()<<p_variable->Name()<<"="<<v
		   <<" vs. {"<<m_bounds[l].first
		   <<","<<m_bounds[l].second<<"}\n";
#endif
    bool res(v<m_bounds[l].first || v>m_bounds[l].second);
    ++l;
    return !m_sel_log->Hit(res);
  }
  if (n==m_nfl[f]) return Trigger(p,l,moms,f+1,0,0);
  moms.push_back(Vec4D());
  for (size_t k(m);k<m_sels[f].size();++k) {
#ifdef DEBUG__Variable_Selector
    msg_Debugging()<<"f = "<<f<<", n = "<<n<<", m = "<<m
		   <<", k = "<<k<<" -> "<<m_cfl[f].IDName()
		   <<" ("<<m_sels[f][k]<<") {\n";
#endif
    moms.back()=m_moms[f][k];
    if (!Trigger(p,l,moms,f,n+1,k+1)) return false;
#ifdef DEBUG__Variable_Selector
    msg_Debugging()<<"}\n";
#endif
  }
  moms.pop_back();
  return true;
}

bool Variable_Selector::Trigger(const Vec4D *p) 
{
#ifdef DEBUG__Variable_Selector
  msg_Debugging()<<METHOD<<"(): {\n";
#endif
  for (size_t j(0);j<m_cfl.size();++j) {
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
    if (p_order!=NULL)
      std::sort(m_moms[j].begin(),m_moms[j].end(),*p_order);
  }
  size_t l(0);
  std::vector<Vec4D> moms;
  bool hit(Trigger(p,l,moms,0,0,0));
#ifdef DEBUG__Variable_Selector
  msg_Debugging()<<"}\n";
#endif
  return hit;
}

