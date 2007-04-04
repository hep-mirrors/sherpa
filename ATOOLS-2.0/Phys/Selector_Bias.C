#include "Selector_Bias.H"

#include "Exception.H"

using namespace ATOOLS;

class Order_Up_E {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return dabs(a[0])>dabs(b[0]); }
};

class Order_Up_ET {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return a.EPerp()>b.EPerp(); }
};

class Order_Up_PT {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return a.PPerp2()>b.PPerp2(); }
};

class Order_Up_Eta {
public:
  bool operator()(const Vec4D &a,const Vec4D &b) const 
  { return dabs(a.Eta())>dabs(b.Eta()); }
};

ET_Bias::ET_Bias(int nin,int nout,Flavour * flavs,ordering::code mode)
{
  m_name = "ET_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_mode = mode;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_sel_log=NULL;
}

ET_Bias::~ET_Bias() 
{
  if (m_fl) delete m_fl;
}

void ET_Bias::BuildCuts(Cut_Data *)
{
}

void ET_Bias::SetRange(std::vector<Flavour> fl,
		       std::vector<std::pair<double,double> > &bd)
{
  if (fl.size()!=1) THROW(fatal_error,"Wrong number of flavours");
  m_bounds=bd;
  m_name="ET_Bias_"+fl.front().IDName();
  m_sels.clear();
  for (int j(m_nin);j<m_n;++j)
    if (fl.front().Includes(m_fl[j])) m_sels.push_back(j);
  m_moms.resize(m_sels.size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool ET_Bias::Trigger(const Vec4D * p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_sels.size();++i) m_moms[i]=p[m_sels[i]]; 
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
  for (size_t i(0);i<Min(m_bounds.size(),m_moms.size());++i) {
    double et(m_moms[i].EPerp());
    msg_Debugging()<<"  "<<i<<" et="<<et<<" vs. {"
		   <<m_bounds[i].first<<","<<m_bounds[i].second<<"}\n";
    if (m_sel_log->Hit(et<m_bounds[i].first ||
		       et>m_bounds[i].second)) return false;
  }
  msg_Debugging()<<"}\n";
  return true;
}

PT_Bias::PT_Bias(int nin,int nout,Flavour * flavs,ordering::code mode)
{
  m_name = "PT_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_mode = mode;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_sel_log=NULL;
}

PT_Bias::~PT_Bias() 
{
  if (m_fl) delete m_fl;
}

void PT_Bias::BuildCuts(Cut_Data *)
{
}

void PT_Bias::SetRange(std::vector<Flavour> fl,
		       std::vector<std::pair<double,double> > &bd)
{
  if (fl.size()!=1) THROW(fatal_error,"Wrong number of flavours");
  m_bounds=bd;
  m_name="PT_Bias_"+fl.front().IDName();
  m_sels.clear();
  for (int j(m_nin);j<m_n;++j)
    if (fl.front().Includes(m_fl[j])) m_sels.push_back(j);
  m_moms.resize(m_sels.size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}


bool PT_Bias::Trigger(const Vec4D * p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_sels.size();++i) m_moms[i]=p[m_sels[i]]; 
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
  for (size_t i(0);i<Min(m_bounds.size(),m_moms.size());++i) {
    double pt(m_moms[i].PPerp());
    msg_Debugging()<<"  "<<i<<" pt="<<pt<<" vs. {"
		   <<m_bounds[i].first<<","<<m_bounds[i].second<<"}\n";
    if (m_sel_log->Hit(pt<m_bounds[i].first ||
		       pt>m_bounds[i].second)) return false;
  }
  msg_Debugging()<<"}\n";
  return true;
}

Eta_Bias::Eta_Bias(int nin,int nout,Flavour * flavs,ordering::code mode)
{
  m_name = "Eta_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_mode = mode;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_sel_log=NULL;
}

Eta_Bias::~Eta_Bias() 
{
  if (m_fl) delete m_fl;
}

void Eta_Bias::BuildCuts(Cut_Data *)
{
}

void Eta_Bias::SetRange(std::vector<Flavour> fl,
			std::vector<std::pair<double,double> > &bd)
{
  if (fl.size()!=1) THROW(fatal_error,"Wrong number of flavours");
  m_bounds=bd;
  m_name="Eta_Bias_"+fl.front().IDName();
  m_sels.clear();
  for (int j(m_nin);j<m_n;++j)
    if (fl.front().Includes(m_fl[j])) m_sels.push_back(j);
  m_moms.resize(m_sels.size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool Eta_Bias::Trigger(const Vec4D * p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_sels.size();++i) {
    m_moms[i]=p[m_sels[i]]; 
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
  for (size_t i(0);i<Min(m_bounds.size(),m_moms.size());++i) {
    double eta(m_moms[i].Eta());
    msg_Debugging()<<"  "<<i<<" eta="<<eta<<" vs. {"
		   <<m_bounds[i].first<<","<<m_bounds[i].second<<"}\n";
    if (m_sel_log->Hit(eta<m_bounds[i].first ||
		       eta>m_bounds[i].second)) return false;
  }
  msg_Debugging()<<"}\n";
  return true;
}

Mass_Bias::Mass_Bias
(int nin,int nout,Flavour * flavs,ordering::code mode)
{
  m_name = "Mass_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_mode = mode;
  m_sel_log=NULL;
}

Mass_Bias::~Mass_Bias() 
{
  if (m_fl) delete m_fl;
}

void Mass_Bias::BuildCuts(Cut_Data *)
{
}

void Mass_Bias::SetRange
(std::vector<Flavour> fl,std::vector<std::pair<double,double> > &bd)
{
  if (fl.size()!=2) THROW(fatal_error,"Wrong number of flavours");
  m_idf=fl[0]==fl[1];
  m_bounds=bd;
  m_name="Mass_Bias_"+fl[0].IDName()+fl[1].IDName();
  m_sels[0].clear();
  m_sels[1].clear();
  for (int i(m_nin);i<m_n;++i) {
    if (fl[0].Includes(m_fl[i])) m_sels[0].push_back(i);
    if (fl[1].Includes(m_fl[i])) m_sels[1].push_back(i);
  }
  m_moms[0].resize(m_sels[0].size());
  m_moms[1].resize(m_sels[1].size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool Mass_Bias::Trigger(const Vec4D * p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t j(0);j<2;++j)
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
  switch (m_mode) {
  case ordering::ET:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_ET());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_ET());
    break;
  case ordering::PT:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_PT());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_PT());
    break;
  case ordering::E:
  default:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_E());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_E());
    break;
  }
  size_t id(0);
  for (size_t j(0);j<m_moms[0].size();++j) {
    for (size_t i(m_idf?j+1:0);i<m_moms[1].size();++i) {
      double m(sqrt((m_moms[0][j]+m_moms[1][i]).Abs2()));
      msg_Debugging()<<"  "<<j<<"&"<<i<<" -> m="<<m
		     <<" vs. {"<<m_bounds[id].first<<","
		     <<m_bounds[id].second<<"}\n";
      if (m_sel_log->Hit(m<m_bounds[id].first || 
			 m>m_bounds[id].second)) return false;
      if (++id>=m_bounds.size()) break;
    }
    if (id>=m_bounds.size()) break;
  }
  msg_Debugging()<<"}\n";
  return true;
}

Delta_Eta_Bias::Delta_Eta_Bias
(int nin,int nout,Flavour * flavs,ordering::code mode)
{
  m_name = "Delta_Eta_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_mode = mode;
  m_sel_log=NULL;
}

Delta_Eta_Bias::~Delta_Eta_Bias() 
{
  if (m_fl) delete m_fl;
}

void Delta_Eta_Bias::BuildCuts(Cut_Data *)
{
}

void Delta_Eta_Bias::SetRange
(std::vector<Flavour> fl,std::vector<std::pair<double,double> > &bd)
{
  if (fl.size()!=2) THROW(fatal_error,"Wrong number of flavours");
  m_idf=fl[0]==fl[1];
  m_bounds=bd;
  m_name="Delta_Eta_Bias_"+fl[0].IDName()+fl[1].IDName();
  m_sels[0].clear();
  m_sels[1].clear();
  for (int i(m_nin);i<m_n;++i) {
    if (fl[0].Includes(m_fl[i])) m_sels[0].push_back(i);
    if (fl[1].Includes(m_fl[i])) m_sels[1].push_back(i);
  }
  m_moms[0].resize(m_sels[0].size());
  m_moms[1].resize(m_sels[1].size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool Delta_Eta_Bias::Trigger(const Vec4D * p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t j(0);j<2;++j)
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
  switch (m_mode) {
  case ordering::ET:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_ET());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_ET());
    break;
  case ordering::PT:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_PT());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_PT());
    break;
  case ordering::E:
  default:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_E());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_E());
    break;
  }
  size_t id(0);
  for (size_t j(0);j<m_moms[0].size();++j) {
    for (size_t i(m_idf?j+1:0);i<m_moms[1].size();++i) {
      double deta(m_moms[0][j].DEta(m_moms[1][i]));
      msg_Debugging()<<"  "<<j<<"&"<<i<<" -> deta="<<deta
		     <<" vs. {"<<m_bounds[id].first<<","
		     <<m_bounds[id].second<<"}\n";
      if (m_sel_log->Hit(deta<m_bounds[id].first || 
			 deta>m_bounds[id].second)) return false;
      if (++id>=m_bounds.size()) break;
    }
    if (id>=m_bounds.size()) break;
  }
  msg_Debugging()<<"}\n";
  return true;
}

Delta_Phi_Bias::Delta_Phi_Bias
(int nin,int nout,Flavour * flavs,ordering::code mode)
{
  m_name = "Delta_Phi_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_mode = mode;
  m_sel_log=NULL;
}

Delta_Phi_Bias::~Delta_Phi_Bias() 
{
  if (m_fl) delete m_fl;
}

void Delta_Phi_Bias::BuildCuts(Cut_Data *)
{
}

void Delta_Phi_Bias::SetRange
(std::vector<Flavour> fl,std::vector<std::pair<double,double> > &bd)
{
  if (fl.size()!=2) THROW(fatal_error,"Wrong number of flavours");
  m_idf=fl[0]==fl[1];
  m_bounds=bd;
  m_name="Delta_Phi_Bias_"+fl[0].IDName()+fl[1].IDName();
  m_sels[0].clear();
  m_sels[1].clear();
  for (int i(m_nin);i<m_n;++i) {
    if (fl[0].Includes(m_fl[i])) m_sels[0].push_back(i);
    if (fl[1].Includes(m_fl[i])) m_sels[1].push_back(i);
  }
  m_moms[0].resize(m_sels[0].size());
  m_moms[1].resize(m_sels[1].size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool Delta_Phi_Bias::Trigger(const Vec4D * p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t j(0);j<2;++j)
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
  switch (m_mode) {
  case ordering::ET:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_ET());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_ET());
    break;
  case ordering::PT:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_PT());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_PT());
    break;
  case ordering::E:
  default:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_E());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_E());
    break;
  }
  size_t id(0);
  for (size_t j(0);j<m_moms[0].size();++j) {
    for (size_t i(m_idf?j+1:0);i<m_moms[1].size();++i) {
      double dphi(m_moms[0][j].DPhi(m_moms[1][i]));
      msg_Debugging()<<"  "<<j<<"&"<<i<<" -> dphi="<<dphi
		     <<" vs. {"<<m_bounds[id].first<<","
		     <<m_bounds[id].second<<"}\n";
      if (m_sel_log->Hit(dphi<m_bounds[id].first || 
			 dphi>m_bounds[id].second)) return false;
      if (++id>=m_bounds.size()) break;
    }
    if (id>=m_bounds.size()) break;
  }
  msg_Debugging()<<"}\n";
  return true;
}

Delta_R_Bias::Delta_R_Bias
(int nin,int nout,Flavour * flavs,ordering::code mode)
{
  m_name = "Delta_R_Bias";
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  m_mode = mode;
  m_sel_log=NULL;
}

Delta_R_Bias::~Delta_R_Bias() 
{
  if (m_fl) delete m_fl;
}

void Delta_R_Bias::BuildCuts(Cut_Data *)
{
}

void Delta_R_Bias::SetRange
(std::vector<Flavour> fl,std::vector<std::pair<double,double> > &bd)
{
  if (fl.size()!=2) THROW(fatal_error,"Wrong number of flavours");
  m_idf=fl[0]==fl[1];
  m_bounds=bd;
  m_name="Delta_R_Bias_"+fl[0].IDName()+fl[1].IDName();
  m_sels[0].clear();
  m_sels[1].clear();
  for (int i(m_nin);i<m_n;++i) {
    if (fl[0].Includes(m_fl[i])) m_sels[0].push_back(i);
    if (fl[1].Includes(m_fl[i])) m_sels[1].push_back(i);
  }
  m_moms[0].resize(m_sels[0].size());
  m_moms[1].resize(m_sels[1].size());
  if (m_sel_log!=NULL) delete m_sel_log;
  m_sel_log = new Selector_Log(m_name);
}

bool Delta_R_Bias::Trigger(const Vec4D * p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t j(0);j<2;++j)
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
  switch (m_mode) {
  case ordering::ET:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_ET());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_ET());
    break;
  case ordering::PT:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_PT());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_PT());
    break;
  case ordering::E:
  default:
    std::sort(m_moms[0].begin(),m_moms[0].end(),Order_Up_E());
    std::sort(m_moms[1].begin(),m_moms[1].end(),Order_Up_E());
    break;
  }
  size_t id(0);
  for (size_t j(0);j<m_moms[0].size();++j) {
    for (size_t i(m_idf?j+1:0);i<m_moms[1].size();++i) {
      double dr(m_moms[0][j].DR(m_moms[1][i]));
      msg_Debugging()<<"  "<<j<<"&"<<i<<" -> dr="<<dr
		     <<" vs. {"<<m_bounds[id].first<<","
		     <<m_bounds[id].second<<"}\n";
      if (m_sel_log->Hit(dr<m_bounds[id].first || 
			 dr>m_bounds[id].second)) return false;
      if (++id>=m_bounds.size()) break;
    }
    if (id>=m_bounds.size()) break;
  }
  msg_Debugging()<<"}\n";
  return true;
}
