#include "Primitive_Integrator.H"

#include "Random.H"
#include "Exception.H"
#include <iomanip>

using namespace ATOOLS;

Primitive_Integrand::~Primitive_Integrand()
{
}

#define DFORMAT std::setw(15)

std::ostream &ATOOLS::operator<<(std::ostream &str,
				 const Primitive_Channel &channel)
{
  str<<"("<<&channel<<"): {\n"
     <<"   m_alpha  = "<<DFORMAT<<channel.Alpha()<<"\n"
     <<"   m_weight = "<<DFORMAT<<channel.Weight()<<"\n"
     <<"   m_sum    = "<<DFORMAT<<channel.Sum()
     <<" -> integral   = "<<DFORMAT<<channel.Mean()<<"\n"
     <<"   m_sum2   = "<<DFORMAT<<channel.Sum2()
     <<" -> error      = "<<DFORMAT<<channel.Sigma()<<"\n"
     <<"   m_np     = "<<DFORMAT<<channel.Points()
     <<" -> rel. error = "<<DFORMAT<<channel.Sigma()/channel.Mean()<<"\n"
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
  m_alpha(0.0), m_weight(0.0), m_sum(0.0), m_sum2(0.0), m_np(0.0) {}

Primitive_Channel::
Primitive_Channel(const std::vector<Primitive_Channel*> &all,
		  Primitive_Channel *const prev,const size_t i):
  m_alpha(0.0), m_weight(0.0), m_sum(0.0), m_sum2(0.0), m_np(0.0)
{
  if (all.empty()) THROW(fatal_error,"No cells.");
  if (i>=prev->m_this.size()) THROW(fatal_error,"Inconsistent dimensions.");
  m_this.resize(prev->m_this.size());
  m_next.resize(prev->m_next.size());
  for (size_t j=0;j<m_this.size();++j) {
    m_this[j]=prev->m_this[j];
    m_next[j]=prev->m_next[j];
  }
  m_this[i]=(m_next[i]->m_this[i]+m_this[i])/2.0;
  m_next[i]=prev->m_next[i];
  prev->m_next[i]=this;
  prev->SetWeight();
  SetWeight();
  m_alpha=(prev->m_alpha/=2.0);
  msg_Debugging()<<"Primitive_Channel::Primitive_Channel("
		 <<prev<<","<<i<<"): {\n"<<*prev<<*this<<"}"<<std::endl;
}
    
Primitive_Channel::~Primitive_Channel()
{
}

double Primitive_Channel::Point(const Primitive_Integrand *function,
				std::vector<double> &point)
{
  if (point.size()!=m_this.size())
    THROW(fatal_error,"Inconsistent dimensions.");
  if (Boundary()) THROW(fatal_error,"Boundary cell selected.");
  for (size_t i=0;i<m_this.size();++i) 
    point[i]=m_this[i]+ran.Get()*(m_next[i]->m_this[i]-m_this[i]);
  double cur=(*function)(point);
  ++m_np;
  m_sum+=cur*m_weight;
  m_sum2+=sqr(cur*m_weight);
  return cur;
}

void Primitive_Channel::Reset()
{
  m_sum=m_sum2=m_np=0.0;
}

void Primitive_Channel::SetWeight()
{
  m_weight=1.0;
  for (size_t i=0;i<m_this.size();++i) 
    m_weight*=m_next[i]->m_this[i]-m_this[i];
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
  root->SetWeight();
}

Primitive_Integrator::Primitive_Integrator():
  m_nopt(1000), m_nmax(1000000), m_nmaxopt(250000), m_error(0.01),
  m_sum(0.0), m_sum2(0.0), m_np(0.0), m_lastdim(0), m_mode(0),
  m_vname("I") {}

Primitive_Integrator::~Primitive_Integrator()
{
  while (!m_channels.empty()) {
    delete m_channels.back();
    m_channels.pop_back();
  }
}

void Primitive_Integrator::SetDimension(const size_t dim)
{
  m_min.resize(dim,0.0);
  m_max.resize(dim,1.0);
}

double Primitive_Integrator::Integrate(const Primitive_Integrand *function)
{
  if (m_min.empty() || m_max.empty())
    THROW(fatal_error,"Zero dimensional integral request.");
  if (m_min.size()!=m_max.size())
    THROW(fatal_error,"Inconsistent dimensions.");
  m_point.resize(m_min.size());
  p_function=function;
  Primitive_Channel::CreateRoot(m_min,m_max,m_channels);
  if (msg.LevelIsDebugging()) {
    msg_Debugging()<<"Primitive_Integrator::Integrate("<<function<<"): {\n";
    for (size_t i=0;i<m_channels.size();++i) msg_Debugging()<<*m_channels[i];
    msg_Debugging()<<"}"<<std::endl;
  }
  double error=1.0;
  while (((long unsigned int)m_np)<Min(m_nmaxopt,m_nmax) && error>m_error) {
    for (long unsigned int n=0;n<m_nopt;++n) Point();
    Update();
    error=dabs(Sigma()/Mean());
    msg_Info()<<om::bold<<m_vname<<om::reset<<" = "<<om::blue<<Mean()
	      <<" "<<m_uname<<om::reset<<" +- ( "<<error*Mean()
	      <<" "<<m_uname<<" = "<<om::red<<error*100.0<<" %"
	      <<om::reset<<" ), n = "<<m_np<<std::endl;
    Split();
  }
  double alpha=1.0/(m_channels.size()-m_point.size());
  for (size_t i=0;i<m_channels.size();++i) 
    if (!m_channels[i]->Boundary()) m_channels[i]->SetAlpha(alpha);
  while (((long unsigned int)m_np)<m_nmax && error>m_error) {
    for (long unsigned int n=0;n<m_nopt;++n) Point();
    Update();
    error=dabs(Sigma()/Mean());
    msg_Info()<<om::bold<<m_vname<<om::reset<<" = "<<om::blue<<Mean()
	      <<" "<<m_uname<<om::reset<<" +- ( "<<error*Mean()
	      <<" "<<m_uname<<" = "<<om::red<<error*100.0<<" %"
	      <<om::reset<<" ), n = "<<m_np<<std::endl;
  }
  return m_sum/m_np;
}

void Primitive_Integrator::Update()
{
  m_sum=m_sum2=m_np=0.0;
  for (size_t i=0;i<m_channels.size();++i) {
    if (!m_channels[i]->Boundary()) {
      m_np+=m_channels[i]->Points();
      m_sum+=m_channels[i]->Mean();
      m_sum2+=m_channels[i]->Variance();
    }
  }
  m_sum2=(m_np-1.0)*m_sum2+m_sum*m_sum*m_np;
  m_sum*=m_np;
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
  selected->Point(p_function,m_point);
}

void Primitive_Integrator::Split()
{
  if (msg.LevelIsDebugging()) {
    msg_Debugging()<<"Primitive_Integrator::Split(): {\n";
    {
      msg_Indent();
      for (size_t i=0;i<m_channels.size();++i) msg_Debugging()<<*m_channels[i];
    }
    msg_Debugging()<<"}"<<std::endl;
  }
  double max=-std::numeric_limits<double>::max(), cur=0.0;
  Primitive_Channel *selected=NULL;
  for (size_t i=0;i<m_channels.size();++i) {
    switch (m_mode) {
    case 1: 
      if ((cur=m_channels[i]->Variance())>max) {
	max=cur;
	selected=m_channels[i];
      }
      break;
    case 0:
    default:
      if ((cur=m_channels[i]->Mean())>max) {
	max=cur;
	selected=m_channels[i];
      }
      break;
    }
  }
  if (selected==NULL) THROW(fatal_error,"Internal error.");
  for (size_t i=0;i<m_channels.size();++i) m_channels[i]->SetAlpha(0.0);
  if (++m_lastdim>=m_point.size()) m_lastdim=0;
  Primitive_Channel *next = 
    new Primitive_Channel(m_channels,selected,m_lastdim);
  m_channels.push_back(next);
  next->SetAlpha(0.5);
  selected->SetAlpha(0.5);
  selected->Reset();
}
