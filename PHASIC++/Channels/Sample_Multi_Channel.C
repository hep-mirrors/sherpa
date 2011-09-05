#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Sample_Multi_Channel.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Channels/VHAAG.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Data_Reader.H"

#include <iomanip>

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {

  std::ostream &operator<<(std::ostream &str,const Idx_Key &key)
  {
    return str<<"{"<<key.m_id<<","<<key.m_t<<","<<key.m_b<<"}";
  }

}

Sample_Multi_Channel::Sample_Multi_Channel
(Phase_Space_Handler *const psh,Color_Integrator *const ci):
  Multi_Channel("SMC_HAAG"), p_colint(ci), m_first(0), m_otfcc(false)
{
  nin=psh->Process()->NIn();
  nout=psh->Process()->NOut();
  if (p_colint==NULL) THROW(fatal_error,"Missing color integrator");
}

Sample_Multi_Channel::~Sample_Multi_Channel()
{
  Clear();
}

void Sample_Multi_Channel::Clear()
{
  while (m_channels.size()) {
    delete m_channels.back();
    m_channels.pop_back();
  }
  while (m_oldchannels.size()) {
    delete m_oldchannels.back();
    m_oldchannels.pop_back();
  }
  m_idx.clear();
}

bool Sample_Multi_Channel::Initialize(const Channel_Vector &chs)
{ 
  int setalpha(2);
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(setalpha,"SMC_SET_ALPHA")) setalpha=2;
  else msg_Info()<<METHOD<<"(): Set setalpha "<<setalpha<<".\n";
  Clear();
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  m_oldchannels=chs;
  VHAAG *firsth(NULL);
  if (m_oldchannels.empty()) {
    if (setalpha&1) setalpha-=1;
  }
  else {
    firsth=(VHAAG*)m_oldchannels.front();
  }
  name = "SMC_HAAG";
  m_mp.resize(nin+nout);
  Idx_Vector perm(nin+nout);
  for (size_t i(0);i<perm.size();++i) perm[i]=i;
  for (int i(1);i<nin+nout;++i) {
    std::vector<size_t> hp(nin+nout);
    for (size_t j(0);j<perm.size();++j) hp[j]=perm[j];
    VHAAG *sc(new VHAAG(nin,nout,hp,firsth));
    if (firsth==NULL) firsth=sc;
    m_channels.push_back(sc);
    if (setalpha&1) {
      double alpha(m_oldchannels[2*Min(i-1,nout-(i-1))]->Alpha());
      sc->SetAlpha(alpha);
    }
    if (i<nin+nout-1) std::swap<Idx_Type>(perm[i],perm[i+1]);
  }
  if (!(setalpha&1))
    for(size_t i(0);i<m_channels.size();++i)
      m_channels[i]->Reset(1.0/m_channels.size());
  Print();
  m_otfcc=p_colint->OTFCC();
  if (!m_otfcc) {
    int otfcc(0);
    if (read.ReadFromFile(otfcc,"CI_OTF_CC")) {
      p_colint->SetOTFCC(m_otfcc=otfcc);
      msg_Info()<<METHOD<<"(): Set on the flight configuration summation "
		<<m_otfcc<<".\n";
    }
  }
  if (m_otfcc==0) FillIdxMap();
  msg_Debugging()<<"} -> "<<name<<"\n";
  if (setalpha>1) {
    Double_Vector alpha(nin+nout-1);
    for (int i(0);i<nin+nout-1;++i) alpha[i]=m_channels[i]->Alpha();
    p_colint->SetAlpha(alpha);
    msg_Tracking()<<"set alpha mode "<<(setalpha>>1)<<"\n";
    p_colint->SetAlphaMode(setalpha>>1);
  }
  return true;
}

void Sample_Multi_Channel::Reset() 
{
  // Multi_Channel::Reset(value);
  Print();
}

Idx_Key Sample_Multi_Channel::GenerateKey(Idx_Vector &perm)
{
  size_t zero(0), one(0);
  for (;zero<perm.size();++zero) if (perm[zero]==0) break;
  Idx_Vector rp(perm.size());
  for (size_t i(0);i<perm.size();++i)
    rp[i]=i+zero<rp.size()?perm[i+zero]:perm[i+zero-rp.size()];
  for (;one<perm.size();++one) if (rp[one]==1) break;
  Idx_Vector mt(rp.size()), mb(rp.size());
  for (size_t i(0);i<2;++i) mt[i]=mb[i]=i;
  for (size_t k(2), i(0);i<rp.size();++i) 
    if (rp[i]>1) {
      mb[mt[k]=rp[i]]=k; 
      ++k;
    }
  for (zero=0;zero<perm.size();++zero) 
    if (perm[zero]==m_first) break;
  for (size_t i(0);i<perm.size();++i)
    rp[i]=i+zero<rp.size()?perm[i+zero]:perm[i+zero-rp.size()];
  perm=rp;
  return Idx_Key(one-1,mt,mb);
}

void Sample_Multi_Channel::FillIdxMap()
{
  while(!p_colint->GeneratePoint());
  m_first=p_colint->Orders().front().front();
  Permutation pp(nin+nout-1);
  Idx_Vector perm(nin+nout);
  for (int j(0);j<pp.MaxNumber();++j) {
    int *pm(pp.Get(j));
    perm[0]=perm.size()-1;
    for (size_t i(1);i<perm.size();++i) perm[i]=pm[i-1];
    Idx_Key key(GenerateKey(perm));
    m_idx[perm]=key;
    msg_Debugging()<<"map "<<perm<<" => "<<m_idx[perm]<<"\n";
  }
}

double Sample_Multi_Channel::SingleWeight
(ATOOLS::Vec4D *p,Cut_Data *cuts,const Idx_Key &key)
{
  msg_Debugging()<<"weight "<<key<<": "<<std::flush;
  for (int j(0);j<nin+nout;++j) m_mp[j]=p[key.m_t[j]];
  m_channels[key.m_id]->GenerateWeight(&m_mp.front(),cuts);
  msg_Debugging()<<m_channels[key.m_id]->Name()<<" "
		 <<m_channels[key.m_id]->Weight()<<"\n";
  return m_channels[key.m_id]->Weight();
}

void Sample_Multi_Channel::OTFWeight
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  m_weight=0.0;
  size_t i(0);
  while (p_colint->NextOrder()) {
    Idx_Vector perm(p_colint->Orders().front());
    if (i++!=m_cur) {
      Idx_Key key(GenerateKey(perm));
      if (m_channels[key.m_id]->Alpha()>0.0)
	m_weight+=m_channels[key.m_id]->Alpha()/
	  SingleWeight(p,cuts,key);
    }
    if (i++!=m_cur) {
      Idx_Vector rperm(perm.size());
      rperm.front()=perm.front();
      for (size_t j(1);j<rperm.size();++j) 
	rperm[j]=perm[perm.size()-j];
      Idx_Key rkey(GenerateKey(rperm));
      if (m_channels[rkey.m_id]->Alpha()>0.0)
	m_weight+=m_channels[rkey.m_id]->Alpha()/
	  SingleWeight(p,cuts,rkey);
    }
  }
  m_weight+=m_channels[m_confs.front().m_id]->Alpha()/
    SingleWeight(p,cuts,m_confs.front());
  if (m_weight!=0.0) m_weight=m_wsum/m_weight;
  msg_Debugging()<<METHOD<<"(): m_weight = "<<m_weight<<"\n";
}

void Sample_Multi_Channel::PDWeight
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  m_weight=0.0;
  for (size_t i(0);i<m_confs.size();++i)
    if (i!=m_cur && m_alpha[i]>0.0)
      m_weight+=m_alpha[i]/SingleWeight(p,cuts,m_confs[i]);
  m_weight+=m_alpha[m_cur]/SingleWeight(p,cuts,m_confs[m_cur]);
  if (m_weight!=0.0) m_weight=1.0/m_weight;
  msg_Debugging()<<METHOD<<"(): m_weight = "<<m_weight<<"\n";
}

void Sample_Multi_Channel::GenerateWeight
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  if (m_otfcc) OTFWeight(p,cuts);
  else PDWeight(p,cuts);
}

void Sample_Multi_Channel::OTFPoint
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  double wsum(0.0);
  std::vector<double> asum(0);
  while (p_colint->NextOrder()) {
    Idx_Vector perm(p_colint->Orders().front());
    size_t type(p_colint->IdentifyType(perm));
    asum.push_back(wsum+=m_channels[type]->Alpha());
    Idx_Vector rperm(perm.size());
    rperm.front()=perm.front();
    for (size_t j(1);j<rperm.size();++j) rperm[j]=perm[perm.size()-j];
    size_t rtype(p_colint->IdentifyType(rperm));
    asum.push_back(wsum+=m_channels[rtype]->Alpha());
  }
  size_t l(0), r(asum.size()-1), i((l+r)/2);
  double disc(ran->Get()*wsum), a(asum[i]);
  while (r-l>1) {
    if (disc<a) r=i;
    else l=i;
    i=(l+r)/2;
    a=asum[i];
  }
  if (disc<asum[l]) r=l;
  for (size_t i(0);i<=r/2;++i) 
    if (!p_colint->NextOrder()) THROW(fatal_error,"Index out of bounds");
  p_colint->RestartOrders();
  Idx_Vector perm(p_colint->Orders().front());
  if (r%2==1) {
    Idx_Vector rperm(perm.size());
    rperm.front()=perm.front();
    for (size_t j(1);j<rperm.size();++j) rperm[j]=perm[perm.size()-j];
    perm=rperm;
  }
  m_confs.resize(1);
  m_confs.front()=GenerateKey(perm);
  m_cur=r;
  m_wsum=wsum;
  msg_Debugging()<<"selected "<<r<<" from l="<<asum[l]<<" < d="
		 <<disc<<" < r="<<asum[r]<<" => "<<perm<<" <-> "
		 <<m_confs.front()<<"\n";
  for (int i(0);i<2;++i) m_mp[i]=p[i];
  m_channels[m_confs.front().m_id]->GeneratePoint(&m_mp.front(),cuts);
  for (int i(2);i<nin+nout;++i) p[i]=m_mp[m_confs.front().m_b[i]];
}

void Sample_Multi_Channel::PDPoint
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  double wsum(0.0);
  m_confs.resize(2*p_colint->Orders().size());
  m_alpha.resize(m_confs.size());
  std::vector<double> asum(m_confs.size());
  for (size_t i(0);i<m_confs.size()/2;++i) {
    size_t id(2*i);
    const Idx_Vector &perm(p_colint->Orders()[i]);
    m_confs[id]=m_idx[perm];
    asum[id]=wsum+=m_alpha[id]=m_channels[m_confs[id].m_id]->Alpha();
    msg_Debugging()<<"select "<<i<<" "<<perm<<" "<<m_confs[id]<<": "
		   <<m_channels[m_confs[id].m_id]->Name()<<" "
		   <<m_alpha[id]<<" -> "<<wsum<<"\n";
    ++id;
    Idx_Vector rperm(perm.size());
    rperm.front()=perm.front();
    for (size_t j(1);j<rperm.size();++j) rperm[j]=perm[perm.size()-j];
    m_confs[id]=m_idx[rperm];
    asum[id]=wsum+=m_alpha[id]=m_channels[m_confs[id].m_id]->Alpha();
    msg_Debugging()<<"select "<<i<<" "<<rperm<<" "<<m_confs[id]<<": "
		   <<m_channels[m_confs[id].m_id]->Name()<<" "
		   <<m_alpha[id]<<" -> "<<wsum<<"\n";
  }
  for (size_t i(0);i<m_alpha.size();++i) m_alpha[i]/=wsum; 
  size_t l(0), r(asum.size()-1), i((l+r)/2);
  double disc(ran->Get()*wsum), a(asum[i]);
  while (r-l>1) {
    if (disc<a) r=i;
    else l=i;
    i=(l+r)/2;
    a=asum[i];
  }
  if (disc<asum[l]) r=l;
  m_cur=r;
  const Idx_Key &selected(m_confs[m_cur]);
  msg_Debugging()<<"selected "<<r<<" from l="<<asum[l]<<" < d="
		 <<disc<<" < r="<<asum[r]<<" => "<<selected<<"\n";
  for (int i(0);i<2;++i) m_mp[i]=p[i];
  m_channels[selected.m_id]->GeneratePoint(&m_mp.front(),cuts);
  for (int i(2);i<nin+nout;++i) p[i]=m_mp[selected.m_b[i]];
}

void Sample_Multi_Channel::GeneratePoint
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  if (m_otfcc) OTFPoint(p,cuts);
  else PDPoint(p,cuts);
}

void Sample_Multi_Channel::AddPoint(double value)
{ 
  m_lastdice=-1;
  Multi_Channel::AddPoint(value);
  if (m_otfcc) m_channels[m_confs.front().m_id]->AddPoint(value);
  else m_channels[m_confs[m_cur].m_id]->AddPoint(value);
}

void Sample_Multi_Channel::Optimize(double error)
{ 
  for (size_t i(0);i<m_channels.size();++i)
    m_channels[i]->Optimize();
  Print();
}

void Sample_Multi_Channel::EndOptimize(double error)
{ 
}

void Sample_Multi_Channel::WriteOut(std::string pid)
{ 
  std::ofstream ofile((pid+"_"+name).c_str());
  ofile.precision(12);
  ofile<<m_channels.size()<<" "<<name<<" "<<n_points<<" "
       <<n_contrib<<" "<<m_optcnt<<std::endl;
  for (size_t i(0);i<m_channels.size();++i) 
    ofile<<m_channels[i]->Name()<<" "<<m_channels[i]->N()<<" "
	 <<m_channels[i]->Alpha()<<std::endl;
  ofile.close();
  for (size_t i(0);i<m_channels.size();++i)
    m_channels[i]->WriteOut(pid);
}
    
bool Sample_Multi_Channel::ReadIn(std::string pid)
{
  std::ifstream ifile((pid+"_"+name).c_str());
  if (ifile.bad()) return false;
  size_t rsize(0);
  std::string rname;
  ifile>>rsize>>rname;
  if (rsize!=m_channels.size() || rname!=name) {
    msg_Error()<<METHOD<<"("<<pid<<"): Inconsistent MC file."<<std::endl;
    return false;
  }
  ifile>>n_points>>n_contrib>>m_optcnt;
  long unsigned int points(0);
  double alpha(0.0);
  for (size_t i(0);i<m_channels.size();++i) {
    ifile>>rname>>points>>alpha;
    if (rname!=m_channels[i]->Name()) {
      msg_Error()<<METHOD<<"("<<pid
		 <<"): Inconsistent channel name."<<std::endl;
      return false;
    }
    m_channels[i]->SetN(points);
    m_channels[i]->SetAlpha(alpha);
  }
  for (size_t i(0);i<m_channels.size();++i)
    m_channels[i]->ReadIn(pid);
  return true;
}

std::string Sample_Multi_Channel::Name()
{
  return name;
}

std::string Sample_Multi_Channel::ChID()
{
  return name;
}

void Sample_Multi_Channel::ISRInfo(int &,double &,double &)
{
}

void Sample_Multi_Channel::Print() 
{
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<"   "<<name<<" {\n";
  for (size_t i(0);i<m_channels.size();++i) 
    msg_Out()<<std::setw(6)<<std::right<<i<<": "<<m_channels[i]->Name()
	     <<": "<<std::setw(15)<<m_channels[i]->Alpha()<<"\n";
  msg_Out()<<"   }\n";
}                 

