#ifndef ATOOLS_Phys_Selector_Bias_H
#define ATOOLS_Phys_Selector_Bias_H

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Phys/Ordering.H"

namespace PHASIC {

  class ET_Bias : public Selector_Base {
    ATOOLS::Order_Base *p_order;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<ATOOLS::Vec4D> m_moms;
  public:
    ET_Bias(int,int,ATOOLS::Flavour *,std::string);
    ~ET_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<ATOOLS::Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const ATOOLS::Vec4D_Vector & );
  };

  class PT_Bias : public Selector_Base {
    ATOOLS::Order_Base *p_order;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<ATOOLS::Vec4D> m_moms;
  public:
    PT_Bias(int,int,ATOOLS::Flavour *,std::string);
    ~PT_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<ATOOLS::Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const ATOOLS::Vec4D_Vector & );
  };

  class Eta_Bias : public Selector_Base {
    ATOOLS::Order_Base *p_order;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<ATOOLS::Vec4D> m_moms;
  public:
    Eta_Bias(int,int,ATOOLS::Flavour *,std::string);
    ~Eta_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<ATOOLS::Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const ATOOLS::Vec4D_Vector & );
  };

  class Mass_Bias : public Selector_Base {
    ATOOLS::Order_Base *p_order;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels[2];
    std::vector<ATOOLS::Vec4D> m_moms[2];
    bool m_idf;
  public:
    Mass_Bias(int,int,ATOOLS::Flavour *,std::string="ET_UP");
    ~Mass_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<ATOOLS::Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const ATOOLS::Vec4D_Vector & );
  };

  class Delta_Eta_Bias : public Selector_Base {
    ATOOLS::Order_Base *p_order;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels[2];
    std::vector<ATOOLS::Vec4D> m_moms[2];
    bool m_idf;
  public:
    Delta_Eta_Bias(int,int,ATOOLS::Flavour *,std::string="ET_UP");
    ~Delta_Eta_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<ATOOLS::Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const ATOOLS::Vec4D_Vector & );
  };

  class Delta_Phi_Bias : public Selector_Base {
    ATOOLS::Order_Base *p_order;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels[2];
    std::vector<ATOOLS::Vec4D> m_moms[2];
    bool m_idf;
  public:
    Delta_Phi_Bias(int,int,ATOOLS::Flavour *,std::string="ET_UP");
    ~Delta_Phi_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<ATOOLS::Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const ATOOLS::Vec4D_Vector & );
  };

  class Delta_R_Bias : public Selector_Base {
    ATOOLS::Order_Base *p_order;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels[2];
    std::vector<ATOOLS::Vec4D> m_moms[2];
    bool m_idf;
  public:
    Delta_R_Bias(int,int,ATOOLS::Flavour *,std::string="ET_UP");
    ~Delta_R_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<ATOOLS::Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const ATOOLS::Vec4D_Vector & );
  };

}// end of namespace ATOOLS

#endif

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;

ET_Bias::ET_Bias(int nin,int nout,Flavour * flavs,std::string mode):
  Selector_Base("ET_Bias")
{
  m_nin  = nin;
  m_nout = nout;
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
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

bool ET_Bias::Trigger(const Vec4D_Vector & p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_sels.size();++i) m_moms[i]=p[m_sels[i]]; 
  std::sort(m_moms.begin(),m_moms.end(),*p_order);
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

namespace PHASIC{

DECLARE_ND_GETTER(ET_Bias_Getter,"ET_Bias",Selector_Base,Selector_Key,true);

Selector_Base *ET_Bias_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  std::string values=key[0][1];
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector flavs(1,flav);
  Data_Reader reader(",",":","!","=");
  reader.SetString(values);
  reader.SetString(values);
  std::vector<std::vector<double> > crits;
  std::vector<std::pair<double,double> > bounds;
  if (!reader.MatrixFromString(crits,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
  for (size_t i(0);i<crits.size();++i) {
    if (crits[i].size()<2) 
      THROW(critical_error,"Invalid Syntax in Selector.dat: '"
	    +values+"'");
    bounds.push_back(std::pair<double,double>
		     (crits[i][0],crits[i][1]));
  }
  ET_Bias *sel = new ET_Bias
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front(),key[0][2]);
  sel->SetRange(flavs,bounds);
  return sel;
}

void ET_Bias_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"transverse energy selector"; 
}

}

PT_Bias::PT_Bias(int nin,int nout,Flavour * flavs,std::string mode):
  Selector_Base("PT_Bias")
{
  m_nin  = nin;
  m_nout = nout;
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
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


bool PT_Bias::Trigger(const Vec4D_Vector & p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_sels.size();++i) m_moms[i]=p[m_sels[i]]; 
  std::sort(m_moms.begin(),m_moms.end(),*p_order);
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

namespace PHASIC{

DECLARE_ND_GETTER(PT_Bias_Getter,"PT_Bias",Selector_Base,Selector_Key,true);

Selector_Base *PT_Bias_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  std::string values=key[0][1];
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector flavs(1,flav);
  Data_Reader reader(",",":","!","=");
  reader.SetString(values);
  reader.SetString(values);
  std::vector<std::vector<double> > crits;
  std::vector<std::pair<double,double> > bounds;
  if (!reader.MatrixFromString(crits,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
  for (size_t i(0);i<crits.size();++i) {
    if (crits[i].size()<2) 
      THROW(critical_error,"Invalid Syntax in Selector.dat: '"
	    +values+"'");
    bounds.push_back(std::pair<double,double>
		     (crits[i][0],crits[i][1]));
  }
  PT_Bias *sel = new PT_Bias
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front(),key[0][2]);
  sel->SetRange(flavs,bounds);
  return sel;
}

void PT_Bias_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"transverse momentum selector"; 
}

}

Eta_Bias::Eta_Bias(int nin,int nout,Flavour * flavs,std::string mode):
  Selector_Base("Eta_Bias")
{
  m_nin  = nin;
  m_nout = nout;
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
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

bool Eta_Bias::Trigger(const Vec4D_Vector & p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_sels.size();++i) {
    m_moms[i]=p[m_sels[i]]; 
  }
  std::sort(m_moms.begin(),m_moms.end(),*p_order);
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

namespace PHASIC{

DECLARE_ND_GETTER(Eta_Bias_Getter,"Eta_Bias",Selector_Base,Selector_Key,true);

Selector_Base *Eta_Bias_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<3) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  std::string values=key[0][1];
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav = flav.Bar();
  Flavour_Vector flavs(1,flav);
  Data_Reader reader(",",":","!","=");
  reader.SetString(values);
  reader.SetString(values);
  std::vector<std::vector<double> > crits;
  std::vector<std::pair<double,double> > bounds;
  if (!reader.MatrixFromString(crits,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
  for (size_t i(0);i<crits.size();++i) {
    if (crits[i].size()<2) 
      THROW(critical_error,"Invalid Syntax in Selector.dat: '"
	    +values+"'");
    bounds.push_back(std::pair<double,double>
		     (crits[i][0],crits[i][1]));
  }
  Eta_Bias *sel = new Eta_Bias
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front(),key[0][2]);
  sel->SetRange(flavs,bounds);
  return sel;
}

void Eta_Bias_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"rapidity selector"; 
}

}

Mass_Bias::Mass_Bias
(int nin,int nout,Flavour * flavs,std::string mode):
  Selector_Base("Mass_Bias")
{
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
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

bool Mass_Bias::Trigger(const Vec4D_Vector & p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t j(0);j<2;++j)
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
  std::sort(m_moms[0].begin(),m_moms[0].end(),*p_order);
  std::sort(m_moms[1].begin(),m_moms[1].end(),*p_order);
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

namespace PHASIC{

DECLARE_ND_GETTER(Mass_Bias_Getter,"Mass_Bias",Selector_Base,Selector_Key,true);

Selector_Base *Mass_Bias_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  std::string values=key[0][2];
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav  = flav.Bar();
  Flavour flav2 = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav2 = flav2.Bar();
  Flavour_Vector flavs(1,flav);
  flavs.push_back(flav2);
  Data_Reader reader(",",":","!","=");
  reader.SetString(values);
  reader.SetString(values);
  std::vector<std::vector<double> > crits;
  std::vector<std::pair<double,double> > bounds;
  if (!reader.MatrixFromString(crits,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
  bounds.clear();
  for (size_t i(0);i<crits.size();++i) {
    if (crits[i].size()<2) 
      THROW(critical_error,"Invalid Syntax in Selector.dat: '"
	    +values+"'");
    bounds.push_back(std::pair<double,double>
		     (crits[i][0],crits[i][1]));
  }
  Mass_Bias *sel = new Mass_Bias
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front(),key[0][3]);
  sel->SetRange(flavs,bounds);
  return sel;
}

void Mass_Bias_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"mass selector"; 
}

}

Delta_Eta_Bias::Delta_Eta_Bias
(int nin,int nout,Flavour * flavs,std::string mode):
  Selector_Base("Delta_Eta_Bias")
{
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
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

bool Delta_Eta_Bias::Trigger(const Vec4D_Vector & p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t j(0);j<2;++j)
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
  std::sort(m_moms[0].begin(),m_moms[0].end(),*p_order);
  std::sort(m_moms[1].begin(),m_moms[1].end(),*p_order);  
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

namespace PHASIC{

DECLARE_ND_GETTER(Delta_Eta_Bias_Getter,"Delta_Eta_Bias",Selector_Base,Selector_Key,true);

Selector_Base *Delta_Eta_Bias_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  std::string values=key[0][2];
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav  = flav.Bar();
  Flavour flav2 = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav2 = flav2.Bar();
  Flavour_Vector flavs(1,flav);
  flavs.push_back(flav2);
  Data_Reader reader(",",":","!","=");
  reader.SetString(values);
  reader.SetString(values);
  std::vector<std::vector<double> > crits;
  std::vector<std::pair<double,double> > bounds;
  if (!reader.MatrixFromString(crits,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
  bounds.clear();
  for (size_t i(0);i<crits.size();++i) {
    if (crits[i].size()<2) 
      THROW(critical_error,"Invalid Syntax in Selector.dat: '"
	    +values+"'");
    bounds.push_back(std::pair<double,double>
		     (crits[i][0],crits[i][1]));
  }
  Delta_Eta_Bias *sel = new Delta_Eta_Bias
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front(),key[0][3]);
  sel->SetRange(flavs,bounds);
  return sel;
}

void Delta_Eta_Bias_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"\\Delta\\eta selector"; 
}

}

Delta_Phi_Bias::Delta_Phi_Bias
(int nin,int nout,Flavour * flavs,std::string mode):
  Selector_Base("Delta_Phi_Bias")
{
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
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

bool Delta_Phi_Bias::Trigger(const Vec4D_Vector & p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t j(0);j<2;++j)
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
  std::sort(m_moms[0].begin(),m_moms[0].end(),*p_order);
  std::sort(m_moms[1].begin(),m_moms[1].end(),*p_order);  
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

namespace PHASIC{

DECLARE_ND_GETTER(Delta_Phi_Bias_Getter,"Delta_Phi_Bias",Selector_Base,Selector_Key,true);

Selector_Base *Delta_Phi_Bias_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  std::string values=key[0][2];
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav  = flav.Bar();
  Flavour flav2 = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav2 = flav2.Bar();
  Flavour_Vector flavs(1,flav);
  flavs.push_back(flav2);
  Data_Reader reader(",",":","!","=");
  reader.SetString(values);
  reader.SetString(values);
  std::vector<std::vector<double> > crits;
  std::vector<std::pair<double,double> > bounds;
  if (!reader.MatrixFromString(crits,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
  bounds.clear();
  for (size_t i(0);i<crits.size();++i) {
    if (crits[i].size()<2) 
      THROW(critical_error,"Invalid Syntax in Selector.dat: '"
	    +values+"'");
    bounds.push_back(std::pair<double,double>
		     (crits[i][0],crits[i][1]));
  }
  Delta_Phi_Bias *sel = new Delta_Phi_Bias
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front(),key[0][3]);
  sel->SetRange(flavs,bounds);
  return sel;
}

void Delta_Phi_Bias_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"\\Delta\\phi selector"; 
}

}

Delta_R_Bias::Delta_R_Bias
(int nin,int nout,Flavour * flavs,std::string mode):
  Selector_Base("Delta_R_Bias")
{
  m_nin  = nin;
  m_nout = nout;
  m_n    = m_nin+m_nout;
  m_fl   = new Flavour[m_n];
  for (int i(0);i<m_n;++i) m_fl[i] = flavs[i];
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
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

bool Delta_R_Bias::Trigger(const Vec4D_Vector & p) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t j(0);j<2;++j)
    for (size_t i(0);i<m_sels[j].size();++i)
      m_moms[j][i]=p[m_sels[j][i]];
  std::sort(m_moms[0].begin(),m_moms[0].end(),*p_order);
  std::sort(m_moms[1].begin(),m_moms[1].end(),*p_order);
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

namespace PHASIC{

DECLARE_ND_GETTER(Delta_R_Bias_Getter,"Delta_R_Bias",Selector_Base,Selector_Key,true);

Selector_Base *Delta_R_Bias_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  int crit1=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][0]));
  int crit2=ToType<int>(key.p_read->Interpreter()->Interprete(key[0][1]));
  std::string values=key[0][2];
  Flavour flav = Flavour((kf_code)abs(crit1));
  if (crit1<0) flav  = flav.Bar();
  Flavour flav2 = Flavour((kf_code)abs(crit2));
  if (crit2<0) flav2 = flav2.Bar();
  Flavour_Vector flavs(1,flav);
  flavs.push_back(flav2);
  Data_Reader reader(",",":","!","=");
  reader.SetString(values);
  reader.SetString(values);
  std::vector<std::vector<double> > crits;
  std::vector<std::pair<double,double> > bounds;
  if (!reader.MatrixFromString(crits,""))
    THROW(critical_error,"Invalid Syntax in Selector.dat: '"+values+"'");
  bounds.clear();
  for (size_t i(0);i<crits.size();++i) {
    if (crits[i].size()<2) 
      THROW(critical_error,"Invalid Syntax in Selector.dat: '"
	    +values+"'");
    bounds.push_back(std::pair<double,double>
		     (crits[i][0],crits[i][1]));
  }
  Delta_R_Bias *sel = new Delta_R_Bias
    (key.p_proc->NIn(),key.p_proc->NOut(),
     (Flavour*)&key.p_proc->Process()->Flavours().front(),key[0][3]);
  sel->SetRange(flavs,bounds);
  return sel;
}

void Delta_R_Bias_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"\\Delta R selector"; 
}

}
