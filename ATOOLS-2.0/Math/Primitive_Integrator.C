#include "Primitive_Integrator.H"

#include "Random.H"
#include "Exception.H"
#include "Run_Parameter.H"
#include <iomanip>
#include <unistd.h>

using namespace ATOOLS;

Primitive_Integrand::~Primitive_Integrand()
{
}

#define DFORMAT std::setw(15)

std::ostream &ATOOLS::operator<<(std::ostream &str,
				 const Primitive_Channel &channel)
{
  str<<"("<<&channel<<"): {\n"
     <<"   m_alpha  = "<<DFORMAT<<channel.Alpha()
     <<" <- oldalpha   = "<<DFORMAT<<channel.OldAlpha()<<"\n"
     <<"   m_weight = "<<DFORMAT<<channel.Weight()<<"\n"
     <<"   m_sum    = "<<DFORMAT<<channel.Sum()
     <<" -> integral   = "<<DFORMAT<<channel.Mean()<<"\n"
     <<"   m_ssum   = "<<DFORMAT<<channel.SSum()
     <<" -> integral   = "<<DFORMAT<<channel.SMean()<<"\n"
     <<"   m_sum2   = "<<DFORMAT<<channel.Sum2()
     <<" -> error      = "<<DFORMAT<<channel.Sigma()<<"\n"
     <<"   m_ssum2  = "<<DFORMAT<<channel.SSum2()
     <<" -> error      = "<<DFORMAT<<channel.SSigma()<<"\n"
     <<"   m_np     = "<<DFORMAT<<channel.Points()
     <<" -> rel. error = "<<DFORMAT<<channel.Sigma()/channel.Mean()<<"\n"
     <<"   m_snp    = "<<DFORMAT<<channel.SPoints()
     <<" -> rel. error = "<<DFORMAT<<channel.SSigma()/channel.SMean()<<"\n"
     <<"   m_this   = ";
  for (size_t i=0;i<channel.m_this.size();++i) 
    str<<DFORMAT<<channel.m_this[i]<<" ";
  str<<"\n   m_next   = ";
  for (size_t i=0;i<channel.m_next.size();++i) 
    str<<DFORMAT<<channel.m_next[i]<<" ";
  str<<"\n              ";
  for (size_t i=0;i<channel.m_next.size();++i) 
    str<<DFORMAT<<(channel.m_next[i]!=NULL?
		   channel.m_next[i]->m_this[i]:0.0)<<" ";
  return str<<"\n}"<<std::endl;
}

Primitive_Channel::Primitive_Channel():
  m_alpha(0.0), m_oldalpha(0.0), m_weight(0.0), 
  m_sum(0.0), m_sum2(0.0), m_max(0.0), m_np(0.0),
  m_ssum(0.0), m_ssum2(0.0), m_snp(0.0), m_pos(0) {}

Primitive_Channel::
Primitive_Channel(Primitive_Channel *const prev,const size_t i,
		  const double &pos):
  m_alpha(0.0), m_oldalpha(0.0), m_weight(0.0), 
  m_sum(0.0), m_sum2(0.0), m_max(0.0), m_np(0.0), 
  m_ssum(0.0), m_ssum2(0.0), m_snp(0.0), m_pos(0)
{
  if (prev->Boundary()) THROW(fatal_error,"Attempt to split boundary cell.");
  if (i>=prev->m_this.size()) THROW(fatal_error,"Inconsistent dimensions.");
  m_this.resize(prev->m_this.size());
  m_next.resize(prev->m_next.size());
  for (size_t j=0;j<m_this.size();++j) {
    m_this[j]=prev->m_this[j];
    m_next[j]=prev->m_next[j];
  }
  if (pos!=std::numeric_limits<double>::max()) m_this[i]=pos;
  else m_this[i]=(m_next[i]->m_this[i]+m_this[i])/2.0;
  m_next[i]=prev->m_next[i];
  prev->m_next[i]=this;
  prev->SetWeight();
  SetWeight();
  m_alpha=(prev->m_alpha/=2.0);
  msg_Debugging()<<"Primitive_Channel::Primitive_Channel("
		 <<prev<<","<<i<<"): {\n"<<*prev<<*this<<"}"<<std::endl;
  if (m_weight<=0.0 || prev->m_weight<=0.0) 
    THROW(fatal_error,"New cell has nonpositive weight.");
}
    
Primitive_Channel::~Primitive_Channel()
{
}

double Primitive_Channel::
Point(const Primitive_Integrand *function,
      std::vector<double> &point,std::vector<double> &opt)
{
  if (opt.size()<m_pos+8*m_this.size())
    THROW(fatal_error,"Inconsistent dimensions.");
  double cur=Point(function,point);
  cur*=m_weight;
  for (size_t i=0;i<m_this.size();++i) {
    double mid=(m_this[i]+m_next[i]->m_this[i])/2.0;
    size_t pos=8*(m_this.size()*m_pos+i)+4*(size_t)(point[i]<mid);
    ++opt[pos];
    opt[pos+1]+=cur;
    opt[pos+2]+=sqr(cur);
    opt[pos+3]=ATOOLS::Max(opt[pos+3],cur);
  }
  return cur;
}

double Primitive_Channel::Point(const Primitive_Integrand *function,
				std::vector<double> &point)
{
  if (point.size()!=m_this.size())
    THROW(fatal_error,"Inconsistent dimensions.");
  if (Boundary()) THROW(fatal_error,"Boundary cell selected.");
  for (size_t i=0;i<m_this.size();++i) 
    point[i]=m_this[i]+ran.Get()*(m_next[i]->m_this[i]-m_this[i]);
  double cur=(*function)(point), weight=cur*m_weight;
  if (!(cur<0.0) && !(cur>=0.0)) THROW(critical_error,"Integrand is nan.");
  ++m_np;
  m_sum+=weight;
  m_sum2+=sqr(weight);
  m_max=ATOOLS::Max(m_max,weight);
  return cur;
}

bool Primitive_Channel::Find(const std::vector<double> &point) const
{
  if (point.size()!=m_this.size()) 
    THROW(fatal_error,"Inconsistent dimesions.");
  for (size_t i=0;i<m_this.size();++i) 
    if (point[i]<m_this[i] || 
	(m_next[i]!=NULL && point[i]>=m_next[i]->m_this[i])) return false;
  return true;
}

void Primitive_Channel::Reset()
{
  m_sum=m_sum2=m_max=m_np=0.0;
}

void Primitive_Channel::ResetAll()
{
  m_sum=m_sum2=m_max=m_np=m_ssum=m_ssum2=m_snp=0.0;
}

void Primitive_Channel::SetWeight()
{
  if (m_next[0]==NULL) {
    m_weight=0.0;
    return;
  }
  m_weight=1.0;
  for (size_t i=0;i<m_this.size();++i) 
    m_weight*=m_next[i]->m_this[i]-m_this[i];
}

void Primitive_Channel::SetAlpha(const double &alpha)
{ 
  if (Boundary()) return;
  m_alpha=alpha;
}

void Primitive_Channel::Store(const double &alpha)
{
  if (alpha==0.0) return;
  m_snp+=m_np;
  m_ssum+=m_sum/alpha;
  m_ssum2+=m_sum2/sqr(alpha);
}

bool Primitive_Channel::
WriteOut(std::fstream *const file,
	 std::map<Primitive_Channel*,size_t> &pmap) const
{
  (*file)<<"[ "<<m_alpha<<" "<<m_oldalpha<<" "<<m_sum<<" "
	 <<m_sum2<<" "<<m_max<<" "<<m_np<<" "<<m_ssum<<" "
	 <<m_ssum2<<" "<<m_snp<<" ( ";
  for (size_t i=0;i<m_this.size();++i) (*file)<<m_this[i]<<" ";
  (*file)<<") ( ";
  for (size_t i=0;i<m_next.size();++i) (*file)<<pmap[m_next[i]]<<" ";
  (*file)<<") ]"<<std::endl;
  return true;
}

bool Primitive_Channel::ReadIn(std::fstream *const file,
			       std::map<size_t,Primitive_Channel*> &pmap)
{
  if (file->eof()) return false;
  std::string dummy;
  (*file)>>dummy>>m_alpha>>m_oldalpha>>m_sum>>m_sum2
	 >>m_max>>m_np>>m_ssum>>m_ssum2>>m_snp>>dummy;
  for (size_t i=0;i<m_this.size();++i) {
    if (file->eof()) return false;
    (*file)>>m_this[i];
  }
  (*file)>>dummy>>dummy;
  for (size_t i=0;i<m_next.size();++i) {
    if (file->eof()) return false;
    size_t pos;
    (*file)>>pos;
    m_next[i]=(Primitive_Channel*)pmap[pos];
  }
  (*file)>>dummy>>dummy;
  if (file->eof()) return false;
  return true;
}

void Primitive_Channel::CreateRoot(const std::vector<double> &min,
				   const std::vector<double> &max,
				   std::vector<Primitive_Channel*> &channels)
{
  if (min.size()!=max.size()) THROW(fatal_error,"Inconsistent dimensions.");
  if (!channels.empty()) THROW(fatal_error,"Initialized channels found.");
  Primitive_Channel *root = new Primitive_Channel();
  channels.push_back(root);
  root->m_this=min;
  root->m_next.resize(max.size());
  for (size_t i=0;i<max.size();++i) {
    Primitive_Channel *next = new Primitive_Channel();
    channels.push_back(next);
    next->m_next=std::vector<Primitive_Channel*>(max.size(),NULL);
    next->m_this=min;
    next->m_this[i]=max[i];
    root->m_next[i]=next;
  }
  root->SetAlpha(1.0);
  root->SaveAlpha();
  root->SetWeight();
}

Primitive_Integrator::Primitive_Integrator():
  m_nopt(10000), m_nmax(1000000), m_error(0.01), m_scale (1.0), 
  m_sum(0.0), m_sum2(0.0), m_max(0.0), m_np(0.0), 
  m_ncells(1000), m_lastdim(0), m_mode(0), 
  m_split(1), m_shuffle(1), m_vname("I") {}

Primitive_Integrator::~Primitive_Integrator()
{
  while (!m_channels.empty()) {
    delete m_channels.back();
    m_channels.pop_back();
  }
}

void Primitive_Integrator::SetDimension(const size_t dim)
{
  m_rmin.resize(dim,0.0);
  m_rmax.resize(dim,1.0);
  m_nosplit.clear();
  for (size_t i=0;i<dim;++i) m_nosplit[i]=false;
}

void Primitive_Integrator::Reset()
{
  m_sum=m_sum2=m_max=m_np=0.0; 
  while (!m_channels.empty()) {
    delete m_channels.back();
    m_channels.pop_back();
  }
}

void Primitive_Integrator::Initialize()
{
  Reset();
  if (m_rmin.empty() || m_rmax.empty())
    THROW(fatal_error,"Zero dimensional integral request.");
  if (m_rmin.size()!=m_rmax.size())
    THROW(fatal_error,"Inconsistent dimensions.");
  m_point.resize(m_rmin.size());
  Primitive_Channel::CreateRoot(m_rmin,m_rmax,m_channels);
}

double Primitive_Integrator::Integrate(const Primitive_Integrand *function)
{
  if (m_channels.empty()) Initialize();
  p_function=function;
  if (msg.LevelIsDebugging()) {
    msg_Debugging()<<"Primitive_Integrator::Integrate("<<function<<"): {\n";
    for (size_t i=0;i<m_channels.size();++i) msg_Debugging()<<*m_channels[i];
    msg_Debugging()<<"}"<<std::endl;
  }
  m_opt.resize(16*m_point.size());
  for (size_t i=0;i<16*m_point.size();++i) m_opt[i]=0.0;
  long unsigned int nfirst=(m_channels.size()-m_point.size())*m_nopt/2;
  for (long unsigned int n=0;n<nfirst;++n) Point();
  Split();
  for (size_t i=0;i<m_channels.size();++i) 
    if (m_channels[i]->Alpha()!=0.5) m_channels[i]->SetAlpha(0.0);
  while (((long unsigned int)m_np)<m_nmax &&
	 m_channels.size()-m_point.size()<m_ncells) {
    for (long unsigned int n=0;n<m_nopt;++n) Point();
    double alpha=1.0/(m_channels.size()-m_point.size());
    for (size_t i=0;i<m_channels.size();++i) m_channels[i]->SetAlpha(alpha);
    CheckTime();
    if (Update(0)<m_error) break;
    Split();
  }
  for (long unsigned int n=0;n<m_nopt;++n) Point();
  m_opt.clear();
  double alpha=1.0/(m_channels.size()-m_point.size());
  for (size_t i=0;i<m_channels.size();++i) 
    if (!m_channels[i]->Boundary()) {
      m_channels[i]->SetAlpha(alpha);
      m_channels[i]->Store(alpha);
      m_channels[i]->Reset();
    }
  size_t add=0;
  while (((long unsigned int)m_np)<2.0*m_nmax) {
    for (long unsigned int n=0;n<m_point.size()*m_nopt;++n) Point();
    CheckTime();
    if (Update(1)<m_error && 
	add++>=ATOOLS::Max((size_t)5,m_point.size())) break;
    Shuffle();
  }
  msg_Info()<<std::endl;
  return m_sum/m_np;
}

void Primitive_Integrator::CheckTime() const 
{
  if (rpa.gen.CheckTime()) return;
  msg.Error()<<om::bold<<"Primitive_Integrator::Integrate(..): "
	     <<om::reset<<om::red<<"Timeout. Interrupt integration."
	     <<om::reset<<std::endl;
  kill(getpid(),SIGINT);
}

double Primitive_Integrator::Update(const int mode)
{
  double sum=0.0;
  m_sum=m_sum2=m_max=m_np=0.0;
  for (size_t i=0;i<m_channels.size();++i) {
    if (!m_channels[i]->Boundary()) {
      double alpha=m_channels[i]->Alpha();
      sum+=alpha;
      switch (mode) {
      case 0:
	if (m_channels[i]->Points()<m_nopt/3)
	  msg.Error()<<"Primitive_Integrator::Update(): "
		     <<"Few points in cell. Increase NOpt."<<std::endl;
	m_np+=m_channels[i]->Points();
	m_sum+=m_channels[i]->Sum()/alpha;
	m_sum2+=m_channels[i]->Sum2()/sqr(alpha);
	break;
      case 1:
	m_np+=m_channels[i]->SPoints();
	m_sum+=m_channels[i]->SSum();
	m_sum2+=m_channels[i]->SSum2();
	break;
      }
      m_max=ATOOLS::Max(m_max,m_channels[i]->Max());
    }
  }
  if (!IsEqual(sum,1.0)) 
    THROW(critical_error,"Summation does not agree.");
  double error=dabs(Sigma()/Mean());
  msg_Info()<<"  "<<om::bold<<m_vname<<om::reset<<" = "<<om::blue
	    <<Mean()*m_scale<<" "<<m_uname<<om::reset<<" +- ( "
	    <<error*Mean()*m_scale<<" "<<m_uname<<" = "<<om::red
	    <<error*100.0<<" %"<<om::reset<<" ), n = "<<m_np<<", "
	    <<m_channels.size()-m_point.size()
	    <<" cells      "<<bm::cr<<std::flush;
  return error;
}

void Primitive_Integrator::Point()
{
  Primitive_Channel *selected=NULL;
  double disc=ran.Get(), sum=0.0;
  for (size_t i=0;i<m_channels.size();++i) {
    sum+=m_channels[i]->Alpha();
    if (sum>=disc) {
      selected=m_channels[i];
      break;
    }
  }
  if (selected==NULL) THROW(fatal_error,"No channel selected.");
  if (m_opt.size()>0) selected->Point(p_function,m_point,m_opt);
  else selected->Point(p_function,m_point);
}

void Primitive_Integrator::Point(std::vector<double> &x)
{
  if (x.size()!=m_point.size()) 
    THROW(fatal_error,"Inconsistent dimensions.");
  Point();
  x=m_point;
}

double Primitive_Integrator::Weight(const std::vector<double> &x) const
{
  for (size_t i=0;i<m_channels.size();++i) {
    if (m_channels[i]->Find(x)) {
      return m_channels[i]->Weight()/m_channels[i]->Alpha();
    }
  }
  msg.Error()<<"Primitive_Integrator::Weight(..): "
	     <<"Point out of range."<<std::endl;
  return 0.0;
}

void Primitive_Integrator::Split()
{
  if (m_split==0) return;
  if (msg.LevelIsDebugging()) {
    msg_Debugging()<<"Primitive_Integrator::Split(): {\n";
    {
      msg_Indent();
      for (size_t i=0;i<m_channels.size();++i) 
	msg_Debugging()<<*m_channels[i];
    }
    msg_Debugging()<<"}"<<std::endl;
  }
  double max=-std::numeric_limits<double>::max(), cur=0.0;
  Primitive_Channel *selected=NULL;
  for (size_t i=0;i<m_channels.size();++i) {
    m_channels[i]->SetAlpha(0.0);
    switch (m_mode) {
    case 1:
      if ((cur=m_channels[i]->Max())>max) {
	max=cur;
	selected=m_channels[i];
      }
      break;
    case 0: 
    default:
      if ((cur=m_channels[i]->Variance())>max) {
	max=cur;
	selected=m_channels[i];
      }
      break;
    }
  }
  if (selected==NULL) THROW(fatal_error,"Internal error.");
  SelectDimension(selected->Position());
  for (size_t i=0;i<16*m_point.size();++i) m_opt[i]=0.0;
  selected->SetAlpha(selected->OldAlpha());
  Primitive_Channel *next = new Primitive_Channel(selected,m_lastdim);
  m_channels.push_back(next);
  next->SaveAlpha();
  next->SetAlpha(0.5);
  next->SetPosition(1);
  selected->SaveAlpha();
  selected->SetAlpha(0.5);
  selected->SetPosition(0);
  selected->ResetAll();
}

void Primitive_Integrator::SelectDimension(const size_t pos)
{
  m_lastdim=0;
  double max=-std::numeric_limits<double>::max(), cur=0.0;
  for (size_t i=0;i<m_point.size();++i) {
    if (m_nosplit[i]==true) continue;
    size_t l=8*(m_point.size()*pos+i);
    switch (m_mode) {
    case 1: 
      cur=ATOOLS::Max(m_opt[l+3],m_opt[l+7]);
      break;
    case 0: 
    default:
      cur=ATOOLS::Max((m_opt[l+2]-m_opt[l+1]*m_opt[l+1]/m_opt[l])
		      /(m_opt[l]-1.0),
		      (m_opt[l+6]-m_opt[l+5]*m_opt[l+5]/m_opt[l+4])
		      /(m_opt[l+4]-1.0));
      break;
    }
    if (cur>max || (IsEqual(cur,max) && ran.Get()>0.5)) {
      max=cur;
      m_lastdim=i;
    }
    if (!(cur<0) && !(cur>=0)) THROW(fatal_error,"Criterion is nan.")
  } 
}

bool Primitive_Integrator::Shuffle()
{
  if (m_shuffle==0) {
    for (size_t i=0;i<m_channels.size();++i) 
      if (!m_channels[i]->Boundary()) 
	m_channels[i]->Store(m_channels[i]->Alpha());
    return true;
  }
  double norm=0.0;
  for (size_t i=0;i<m_channels.size();++i) {
    if (!m_channels[i]->Boundary()) {
      double alpha=m_channels[i]->Alpha();
      m_channels[i]->Store(alpha);
#ifdef USING__Old_Shufflemode
      alpha=(m_channels[i]->Sum2()>Accu()?
	     sqrt(m_channels[i]->Sum2()/m_channels[i]->Points())/alpha:
	     sqrt(m_channels[i]->SSum2()/m_channels[i]->SPoints()));
#else
      alpha=sqrt(m_channels[i]->SSum2()/m_channels[i]->SPoints());
#endif
      if (!(alpha>=0.0)) THROW(fatal_error,"Invalid weight.");
      m_channels[i]->SetAlpha(alpha);
      m_channels[i]->Reset();
      norm+=alpha;
    }
  }
  for (size_t i=0;i<m_channels.size();++i) {
    if (!m_channels[i]->Boundary()) 
      m_channels[i]->SetAlpha(m_channels[i]->Alpha()/norm);
  }
  return true;
}

bool Primitive_Integrator::WriteOut(const std::string &filename) const
{
  std::fstream *file = new std::fstream(filename.c_str(),std::ios::out);
  if (file->bad()) {
    delete file;
    return false;
  }
  bool result=true;
  file->precision(14);
  (*file)<<m_nopt<<" "<<m_nmax<<" "<<m_ncells<<"\n";
  (*file)<<m_error<<" "<<m_scale<<"\n";
  (*file)<<m_sum<<" "<<m_sum2<<" "<<m_max<<" "<<m_np<<"\n";
  (*file)<<m_rmin.size()<<" ";
  for (size_t i=0;i<m_rmin.size();++i) (*file)<<m_rmin[i]<<" ";
  (*file)<<"\n"<<m_rmax.size()<<" ";
  for (size_t i=0;i<m_rmax.size();++i) (*file)<<m_rmax[i]<<" ";
  (*file)<<"\n"<<m_channels.size()<<" {\n";
  std::map<Primitive_Channel*,size_t> pmap;
  for (size_t i=0;i<m_channels.size();++i) pmap[m_channels[i]]=i+1;
  for (size_t i=0;i<m_channels.size();++i) 
    if (!m_channels[i]->WriteOut(file,pmap)) result=false;
  (*file)<<"}\n"<<m_reserved.size()<<" {\n";
  for (Position_Map::const_iterator rit=m_reserved.begin();
       rit!=m_reserved.end();++rit) 
    (*file)<<rit->first<<" "<<rit->second.first<<" "
	   <<rit->second.second.first<<" "
	   <<rit->second.second.second<<"\n";
  (*file)<<"}\n"<<std::endl;
  delete file;
  return result;
}

bool Primitive_Integrator::ReadIn(const std::string &filename)
{
  msg_Debugging()<<"Primitive_Integrator::ReadIn(\""
		 <<filename<<"\"):"<<std::endl;
  if (!m_channels.empty()) return false;
  std::fstream *file = new std::fstream(filename.c_str(),std::ios::in);
  if (!file->good()) {
    msg_Info()<<"Primitive_Integrator::ReadIn(\""<<filename<<"\"): "
	      <<"Cannot find file."<<std::endl;
    delete file;
    return false;
  }
  std::string dummy;
  file->precision(14);
  (*file)>>m_nopt>>m_nmax>>m_ncells;
  (*file)>>m_error>>m_scale;
  (*file)>>m_sum>>m_sum2>>m_max>>m_np;
  if (file->eof()) {
    delete file;
    return false;
  }
  size_t size;
  (*file)>>size;
  m_rmin.resize(size);
  for (size_t i=0;i<m_rmin.size();++i) (*file)>>m_rmin[i];
  (*file)>>size;
  m_rmax.resize(size);
  for (size_t i=0;i<m_rmax.size();++i) (*file)>>m_rmax[i];
  (*file)>>size>>dummy;
  if (file->eof() || dummy!="{") {
    msg.Error()<<"Primitive_Integrator::ReadIn("<<filename<<"): "
	       <<"Data error.";
    delete file;
    return false;
  }
  Primitive_Channel::CreateRoot(m_rmin,m_rmax,m_channels);
  std::map<size_t,Primitive_Channel*> pmap;
  pmap[0]=NULL;
  for (size_t i=0;i<size;++i) {
    if (i<m_rmin.size()) {
      pmap[i+1]=m_channels[i];
    }
    else {
      if (i<size-1) m_channels.
	push_back(new Primitive_Channel(m_channels[0],m_lastdim));
      pmap[i+1]=m_channels[i];
    }
  }
  bool result=true;
  for (size_t i=0;i<m_channels.size();++i) {
    if (!m_channels[i]->ReadIn(file,pmap)) result=false;
    msg_Tracking()<<*m_channels[i]<<std::endl;
  }
  for (size_t i=0;i<m_channels.size();++i) m_channels[i]->SetWeight();
  (*file)>>dummy;
  if (file->eof() || dummy!="}") {
    msg.Error()<<"Primitive_Integrator::ReadIn("<<filename<<"): "
	       <<"Data error.";
    delete file;
    return false;
  }
  (*file)>>size>>dummy;
  for (size_t i=0;i<size;++i) {
    std::string tag;
    size_t pos, next, ext;
    (*file)>>tag>>pos>>next>>ext;
    m_reserved[tag]=Position_Pair(pos,std::pair<size_t,size_t>(next,ext));
  }
  (*file)>>dummy;
  if (file->eof() || dummy!="}") {
    msg.Error()<<"Primitive_Integrator::ReadIn("<<filename<<"): "
	       <<"Data error.";
    delete file;
    return false;
  }
  delete file;
  m_point.resize(m_rmin.size());
  return result;
}

void Primitive_Integrator::Reserve(const std::string &key,const size_t n,
				   const size_t nprev)
{
  if (m_reserved.find(key)!=m_reserved.end())
    THROW(critical_error,"Key already present.");
  size_t cur=0;
  for (Position_Map::iterator rit=m_reserved.begin();
       rit!=m_reserved.end();++rit) cur+=rit->second.second.second;
  if (cur+n>m_rmin.size()) THROW(fatal_error,"Inconsistent dimesions.");
  m_reserved[key]=Position_Pair(cur,std::pair<size_t,size_t>(nprev,n));
  for (size_t i=cur;i<cur+nprev;++i) m_nosplit[i]=true;
}

const double *const Primitive_Integrator::
Reserved(const std::string &key) const
{
  Position_Map::const_iterator rit=m_reserved.find(key);
  if (rit==m_reserved.end()) THROW(critical_error,"Key not found.");
  return &m_point[rit->second.first];
}

void Primitive_Integrator::
Split(const std::string &key,
      const size_t nprev,const std::vector<double> &pos)
{
  if (m_channels.empty()) 
    THROW(fatal_error,"No cells. Call Initialize() first.");
  Position_Map::const_iterator rit=m_reserved.find(key);
  if (rit==m_reserved.end()) THROW(critical_error,"Key not found.");
  if (nprev>=rit->second.second.first)
    THROW(critical_error,"Inconsistent dimension.");
  size_t dim=rit->second.first+nprev;
  const std::vector<Primitive_Channel*> channels(m_channels);
  for (size_t i=pos.size();i>0;--i) {
    for (size_t j=0;j<channels.size();++j) {
      if (channels[j]->Boundary()) continue;
      Primitive_Channel *next = 
	new Primitive_Channel(channels[j],dim,pos[i-1]);
      m_channels.push_back(next);    
    }
  }
  const double alpha=1.0/(m_channels.size()-m_point.size());
  for (size_t i=0;i<m_channels.size();++i) {
    m_channels[i]->SetAlpha(alpha);
    m_channels[i]->SaveAlpha();
  }
}

void Primitive_Integrator::SetMin(const std::vector<double> &min) 
{ 
  m_rmin=min; 
  m_nosplit.clear();
  for (size_t i=0;i<min.size();++i) m_nosplit[i]=false;
}

void Primitive_Integrator::SetMax(const std::vector<double> &max) 
{
  m_rmax=max; 
  m_nosplit.clear();
  for (size_t i=0;i<max.size();++i) m_nosplit[i]=false;
}

