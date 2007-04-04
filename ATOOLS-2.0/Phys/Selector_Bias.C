#include "Selector_Bias.H"

using namespace ATOOLS;

class Order_Up_E {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return dabs(a[0])>dabs(b[0]); }
  bool operator()(const std::pair<Vec4D,Vec4D> &a,const std::pair<Vec4D,Vec4D> &b) const 
  { 
    double epa1(a.first[0]), epb1(b.first[0]);
    if (epa1>epb1) return true;
    if (epa1<epb1) return false;
    return a.second[0]>b.second[0];
  }
};

class Order_Up_ET {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return a.EPerp()>b.EPerp(); }
  bool operator()(const std::pair<Vec4D,Vec4D> &a,const std::pair<Vec4D,Vec4D> &b) const 
  { 
    double epa1(a.first.EPerp()), epb1(b.first.EPerp());
    if (epa1>epb1) return true;
    if (epa1<epb1) return false;
    return a.second.EPerp()>b.second.EPerp();
  }
};

class Order_Up_PT {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return a.PPerp2()>b.PPerp2(); }
  bool operator()(const std::pair<Vec4D,Vec4D> &a,const std::pair<Vec4D,Vec4D> &b) const 
  { 
    double epa1(a.first.PPerp2()), epb1(b.first.PPerp2());
    if (epa1>epb1) return true;
    if (epa1<epb1) return false;
    return a.second.PPerp2()>b.second.PPerp2();
  }
};

class Order_Up_Eta {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return dabs(a.Eta())>dabs(b.Eta()); }
};


//------------------------------------------------------------------

ET_Bias::ET_Bias(int nin,int nout,Flavour * flavs)
{
  m_name = "ET_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_sel_log=NULL;
}

ET_Bias::~ET_Bias() {
  if (m_fl) delete m_fl;
}

void ET_Bias::SetRange(std::vector<Flavour> fl,
		       std::vector<std::pair<double,double> > &bd)
{
  m_bounds=bd;
  m_name="ET_Bias_";
  m_sels.clear();
  std::set<int> found;
  for (int j(m_nin);j<m_n;++j) {
    for (size_t i(0);i<fl.size();++i) {
      if (found.find(i)!=found.end()) continue;
      if (fl[i].Includes(m_fl[j])) {
	m_sels.push_back(j);
	m_name+=fl[i];
	found.insert(i);
	break;
      }
    }
  }
  m_moms.resize(m_sels.size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool ET_Bias::Trigger(const Vec4D * p) {
  msg.Out()<<"---------------------------------------------------------------------"<<std::endl;
  for (size_t i(0);i<m_sels.size();++i) {
    m_moms[i]=p[m_sels[i]]; 
    msg.Out()<<METHOD<<":"<<m_moms[i]<<" : "<<m_moms[i].EPerp()<<"/"<<m_moms[i].Eta()<<std::endl;
  }
  std::sort(m_moms.begin(),m_moms.end(),Order_Up_ET());
  for (size_t i(0);i<m_bounds.size();++i) {
    double et(m_moms[i].EPerp());
    if (m_sel_log->Hit(et<m_bounds[i].first ||
		       et>m_bounds[i].second)) {
      msg.Out()<<METHOD<<" mom["<<i<<"] = "<<m_moms[i]<<" failed. "<<std::endl;
      return false;
    }
  }
  return true;
}




//---------------------------------------------------------------------

PT_Bias::PT_Bias(int nin,int nout,Flavour * flavs)
{
  m_name = "PT_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_sel_log=NULL;
}

PT_Bias::~PT_Bias() {
  if (m_fl) delete m_fl;
}

void PT_Bias::SetRange(std::vector<Flavour> fl,
		       std::vector<std::pair<double,double> > &bd)
{
  m_bounds=bd;
  m_name="PT_Bias_";
  m_sels.clear();
  std::set<int> found;
  for (int j(m_nin);j<m_n;++j) {
    for (size_t i(0);i<fl.size();++i) {
      if (found.find(i)!=found.end()) continue;
      if (fl[i].Includes(m_fl[j])) {
	m_sels.push_back(j);
	m_name+=fl[i];
	found.insert(i);
	break;
      }
    }
  }
  m_moms.resize(m_sels.size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}


bool PT_Bias::Trigger(const Vec4D * p) {
  for (size_t i(0);i<m_sels.size();++i) {
    m_moms[i]=p[m_sels[i]]; 
  }
  std::sort(m_moms.begin(),m_moms.end(),Order_Up_PT());
  for (size_t i(0);i<m_bounds.size();++i) {
    double pt(m_moms[i].PPerp());
    if (m_sel_log->Hit(pt<m_bounds[i].first ||
		       pt>m_bounds[i].second)) return false;
  }
  return true;
}




//---------------------------------------------------------------------

Eta_Bias::Eta_Bias(int nin,int nout,Flavour * flavs)
{
  m_name = "Eta_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_sel_log=NULL;
}

Eta_Bias::~Eta_Bias() {
  if (m_fl) delete m_fl;
}

void Eta_Bias::SetRange(std::vector<Flavour> fl,
			std::vector<std::pair<double,double> > &bd)
{
  m_bounds=bd;
  m_name="Eta_Bias_";
  m_sels.clear();
  std::set<int> found;
  for (int j(m_nin);j<m_n;++j) {
    for (size_t i(0);i<fl.size();++i) {
      if (found.find(i)!=found.end()) continue;
      if (fl[i].Includes(m_fl[j])) {
	m_sels.push_back(j);
	m_name+=fl[i];
	found.insert(i);
	break;
      }
    }
  }
  m_moms.resize(m_sels.size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool Eta_Bias::Trigger(const Vec4D * p) {
  for (size_t i(0);i<m_sels.size();++i) {
    m_moms[i]=p[m_sels[i]]; 
  }
  std::sort(m_moms.begin(),m_moms.end(),Order_Up_ET());
  for (size_t i(0);i<m_bounds.size();++i) {
    double eta(m_moms[i].Eta());
    if (m_sel_log->Hit(eta<m_bounds[i].first ||
		       eta>m_bounds[i].second)) {
      msg.Out()<<METHOD<<" mom["<<i<<"] = "<<m_moms[i]<<" failed. "<<std::endl;
      return false;
    }
  }
  return true;
}

//---------------------------------------------------------------------

Leading_Delta_Eta_Bias::Leading_Delta_Eta_Bias(int nin,int nout,Flavour * flavs,
					       ordering::code mode)
{
  m_name = "Leading_Delta_Eta_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_mode = mode;
  m_sel_log=NULL;
}

Leading_Delta_Eta_Bias::~Leading_Delta_Eta_Bias() {
  if (m_fl) delete m_fl;
}

void Leading_Delta_Eta_Bias::SetRange(std::vector<Flavour> flpair,
				      std::vector<std::pair<double,double> > &bd)
{
  if (flpair.size()!=2) {
    msg.Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   Wrong number of flavours to be paired : "<<flpair.size()<<"."<<std::endl;
    abort();
  }
  m_bounds=bd;
  m_name="Leading_Delta_Eta_Bias_";
  m_sels.clear();
  std::set<int> found;
  for (int i(m_nin);i<m_n;++i) {
    for (int j(i+1);j<m_n;++j) {
      if ((flpair[0].Includes(m_fl[i]) && flpair[1].Includes(m_fl[j])) ||
	  (flpair[0].Includes(m_fl[j]) && flpair[1].Includes(m_fl[i]))) {
	m_sels.push_back(std::pair<int,int>(i,j));
	m_sels.push_back(std::pair<int,int>(j,i));
	m_name+=m_fl[i]+m_fl[j];
      }
    }
  }
  m_moms.resize(m_sels.size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool Leading_Delta_Eta_Bias::Trigger(const Vec4D * p) {
  for (size_t i(0);i<m_sels.size();++i) {
    msg_Debugging()<<"fill "<<m_sels[i].first<<","<<m_sels[i].second
		   <<" -> "<<p[m_sels[i].first]<<" "<<p[m_sels[i].second]<<"\n";
    m_moms[i]=std::pair<Vec4D,Vec4D>(p[m_sels[i].first],p[m_sels[i].second]);
  }
  switch (m_mode) {
    case ordering::ET:
      std::sort(m_moms.begin(),m_moms.end(),Order_Up_ET());
      break;
    case ordering::PT:
      std::sort(m_moms.begin(),m_moms.end(),Order_Up_PT());
      break;
    case ordering::E:
    default:
      std::sort(m_moms.begin(),m_moms.end(),Order_Up_E());
      break;
  }
  for (size_t i(0);i<m_bounds.size();++i) {
    double deta(dabs(m_moms[i].first.DEta(m_moms[i].second)));
    if (m_sel_log->Hit(deta<m_bounds[i].first || deta>m_bounds[i].second)) {
      msg.Out()<<METHOD<<" mom["<<i<<"] = {"
	       <<m_moms[i].first<<","<<m_moms[i].second<<"} failed: "<<deta<<std::endl;
      return false;
    }
    else {
      msg.Out()<<METHOD<<" mom["<<i<<"] = {"
	       <<m_moms[i].first<<","<<m_moms[i].second<<"} THROUGH: "<<deta<<std::endl;
    }
  }
  return true;
}

