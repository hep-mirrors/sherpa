#include "Parni.H"
#include "Random.H"
#include <iostream>

#include <vector>
#include <stdlib.h>

using namespace ATOOLS;

template <int t_dim> N_Tuple<t_dim>::N_Tuple() 
{
  for (int i=0;i<t_dim;++i) m_x[i]=0.;
}

template <int t_dim>
N_Tuple<t_dim>::N_Tuple(const N_Tuple<t_dim>& ntup)
{
  for (int i=0;i<t_dim;++i) m_x[i]=ntup[i];  
}

template <int t_dim>
N_Tuple<t_dim> & N_Tuple<t_dim>::operator=(const N_Tuple<t_dim> & ntup)
{
  for (int i=0;i<t_dim;++i) m_x[i]=ntup[i];  
  return *this;
}


template <int t_dim>
bool ATOOLS::operator<(const N_Tuple<t_dim> & a, const N_Tuple<t_dim> & b)
{
  for (int i=0;i<t_dim;++i) if (a[i]>=b[i]) return false;
  return true;
}

template <int t_dim>
bool ATOOLS::operator<=(const N_Tuple<t_dim> & a, const N_Tuple<t_dim> & b)
{
  for (int i=0;i<t_dim;++i) if (a[i]>b[i]) return false;
  return true;
}

template <int t_dim>
bool ATOOLS::operator==(const N_Tuple<t_dim> & a, const N_Tuple<t_dim> & b)
{
  for (int i=0;i<t_dim;++i) if (a[i]!=b[i]) return false;
  return true;
}




template <int t_dim>
Parni_Channel<t_dim>::Parni_Channel() : 
  m_points(0), m_alpha(1.), m_sum(0.), m_sum2(0.), m_min(0.), m_max(0.)
{
  for (int i=0;i<t_dim;++i) m_lower[i]=0.;
  for (int i=0;i<t_dim;++i) m_upper[i]=1.;
}

template <int t_dim>
Parni_Channel<t_dim>::Parni_Channel(const Parni_Channel<t_dim> & pc) :
  m_points(pc.m_points), m_alpha(pc.m_alpha), m_sum(pc.m_sum), m_sum2(pc.m_sum2), m_min(pc.m_min), m_max(pc.m_max)
{
  for (int i=0;i<t_dim;++i) m_lower[i]=pc.m_lower[i];
  for (int i=0;i<t_dim;++i) m_upper[i]=pc.m_upper[i];
}

template <int t_dim>
Parni_Channel<t_dim> & Parni_Channel<t_dim>::operator=(const Parni_Channel<t_dim> & pc)
{
  m_points = pc.m_points; 
  m_alpha  = pc.m_alpha; 
  m_sum  = pc.m_sum; 
  m_sum2 = pc.m_sum2; 
  m_min = pc.m_min;
  m_max = pc.m_max;
  for (int i=0;i<t_dim;++i) m_lower[i]=pc.m_lower[i];
  for (int i=0;i<t_dim;++i) m_upper[i]=pc.m_upper[i];
}

template <int t_dim>
bool Parni_Channel<t_dim>::Includes(const N_Tuple<t_dim> & xy)
{
  return m_include=(m_lower<xy && xy<m_upper);
}

template <int t_dim>
N_Tuple<t_dim> & Parni_Channel<t_dim>::GeneratePoint() 
{
  for (int i=0;i<t_dim;++i) {
    m_point[i]=m_lower[i]+(m_upper[i]-m_lower[i])*ran.Get();
  }
  return m_point;
}

template <int t_dim>
void Parni_Channel<t_dim>::AddPoint(double weight,double alphasum)
{
  ++m_points;

  if (m_min<weight*m_alpha || m_min==-1.) m_min=weight*m_alpha;
  if (m_max<weight*m_alpha) m_max=weight*m_alpha;

//   m_sum+=m_alpha/alphasum*weight;
//   m_sum2+=m_alpha/alphasum*weight*weight;

  m_sum+=weight;
  m_sum2+=weight*weight;
}

template <int t_dim>
double Parni_Channel<t_dim>::OptimizeAlpha()
{
  //  std::cout<<" Parni_Channel<t_dim>::OptimizeAlpha() mode="<<m_mode<<std::endl;
  // mode==8
  m_alpha=m_sum;
  return m_alpha;

  // mode==0
  m_alpha=m_max*Integral();



  //  return m_max*(m_max-m_min)*Integral();
  return m_max*Integral();

  // mode==1  m_alpha*=sqrt(m_sum2/double(m_points));
  Reset();
  return m_alpha;
}

template <int t_dim>
double Parni_Channel<t_dim>::Integral() const
{
  double sum=1;
  for (int i=0; i<t_dim;++i) {
    sum*=(m_upper[i]-m_lower[i]);
  }
  return sum;
}

template <int t_dim>
void Parni_Channel<t_dim>::Reset()
{
  m_points = 0;
  m_sum    = 0.;
  m_sum2   = 0.;
}

template <int t_dim>
bool ATOOLS::operator==(const Parni_Channel<t_dim> & a, const Parni_Channel<t_dim> & b)
{
  return  (a.Upper()==b.Upper() && a.Lower()==b.Lower());
}

template <int t_dim>
std::ostream & ATOOLS::operator<<(std::ostream & s, const Parni_Channel<t_dim> & pc) {
  s<<"Box "<<pc.Alpha()<<" ";
  N_Tuple<t_dim> low=pc.Lower();
  N_Tuple<t_dim> up =pc.Upper();
  s<<"("<<low[0]<<","<<up[0]<<")";
  for (int i=1;i<t_dim;++i) {
    s<<"x("<<low[i]<<","<<up[i]<<")";
  }
  s<<std::endl;
  return s;
}

// ----------------------------------------

template <int t_dim>
Parni<t_dim>::Parni(const std::string & name,int mode)
{
  m_nevt=0;
  m_name=name;
  m_mode=mode;

  if (mode!=8) {
    std::cout<<" Initializing Parni<"<<t_dim<<"> "<<m_name<<std::endl;
    std::cout<<"   create 1 Channels "<<std::endl;
    m_channels.push_back(Parni_Channel<t_dim>());
    SplittChannel(m_channels.begin());
    
    if ((m_mode&1)==0) {
      std::vector<Parni_Channel<t_dim> > store;
      for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) 
	store.push_back(*it);
      m_channels.clear();
      for (int i=0;i<store.size();++i) SplittChannel(store[i],0);
    }
  }
  else {
    SplittChannel();
  }

  double test_sum_alpha=0.;
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    it->SetAlpha(it->Alpha()/2.);
    test_sum_alpha+=it->Alpha();
  }
  std::cout<<" test sum alpha="<<test_sum_alpha<<std::endl;  

  if (!IsEqual(test_sum_alpha,1.)) {
    std::cout<<" norming channels"<<test_sum_alpha<<std::endl;  
    for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
      it->SetAlpha(it->Alpha()/test_sum_alpha);
    }
  }
}

template <int t_dim>
N_Tuple<t_dim> & Parni<t_dim>::GeneratePoint()
{
  double rn=ran.Get();
  double asum=0;
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    asum+=it->Alpha();
    if (rn<=asum) {
      m_selected=it;
      return it->GeneratePoint();
    }
  }
  std::cout<<" Warning in Parni<"<<t_dim<<">  rn="<<rn<<" vs. asum="<<asum<<std::endl;
  std::cout<<"  try again ... "<<std::endl;
  return GeneratePoint();
}

template <int t_dim>
double Parni<t_dim>::GenerateWeight(const N_Tuple<t_dim> & xy)
{
  int hit=0;
  m_weight=0.;
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    if (it->Includes(xy)) { 
      //      std::cout<<" Weight+="<<it->Alpha()<<std::endl;
      m_weight+=it->Alpha()/it->Integral();
      ++hit;
    }
  }
  //  if (hit!=1) std::cout<<" found "<<xy[0]<<" int "<<hit<<" boxes "<<std::endl;
  return m_weight;
}

template <int t_dim>
void Parni<t_dim>::AddPoint(double value)
{
  ++m_nevt;
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    if (it->Includes()) { 
      it->AddPoint(value,m_weight);
    }
  }
}

template <int t_dim>
void Parni<t_dim>::AddPoint(const N_Tuple<t_dim> & xy, double value)
{
  ++m_nevt;
  //  std::cout<<" AddPoint("<<xy[0]<<","<<xy[1]<<", w="<<value<<")"<<std::endl;

  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    if (it->Includes(xy)) { 
      it->AddPoint(value,1.);
    }
  }

  if (m_nevt%(t_dim*2500)==0) { 
    double sum_alpha=0.;
    for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
      double alpha=it->OptimizeAlpha();
      sum_alpha+=alpha;
    }
    std::cout<<sum_alpha/m_nevt<<std::endl;
    WriteOut("histo.out");
  }
}


template <int t_dim>
void Parni<t_dim>::Optimize()
{
  std::cout<<" Parni<"<<t_dim<<">::Optimize() with "<<m_channels.size()<<" channels "<<std::endl;
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    std::cout<<" npoints="<<it->NPoints()<<std::endl;
    std::cout<<*it;
  }
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    if ((m_mode&1)==1 && it->NPoints()<=10) return;
    if (it->NPoints()<1) return;
  }

  double sum_alpha=0;
  Channel_Iterator win=m_channels.begin();
  double max_alpha=0.;
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    std::cout<<" "<<it->Alpha()<<" -> ";
    double alpha=it->OptimizeAlpha();
    std::cout<<alpha<<std::endl;
    sum_alpha+=alpha;
    if (alpha>max_alpha) {
      max_alpha=alpha;
      win=it;
    }
  }


  if (m_channels.size()<400) {
    SplittChannel(win);

    sum_alpha+=max_alpha;
  }
  //  if (m_channels.size()>5) Merge();

  sum_alpha=0.;
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    sum_alpha+=it->Alpha();
  }


  double test_sum_alpha=0.;
  double min_alpha=1.;
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    it->SetAlpha(it->Alpha()/sum_alpha);
    if (it->Alpha()<min_alpha) min_alpha=it->Alpha();
    std::cout<<" "<<min_alpha;
    test_sum_alpha+=it->Alpha();
  }
  if (min_alpha*m_channels.size()<1.e-3) {
    std::cout<<" MMMMM "<<std::endl;
    Merge(2);
  }
  std::cout<<" test sum alpha="<<test_sum_alpha<<std::endl;
}

template <int t_dim>
void Parni<t_dim>::SplittChannel()
{
  m_channels.clear();

  N_Tuple<t_dim> ns;
  N_Tuple<t_dim> is;
  int  n=32;
  for (int i=0; i<t_dim; ++i) {
    ns[i]=n;
    is[i]=0;
    n/=2;
  }

  bool complete=false;
  for (;!complete;) {
    N_Tuple<t_dim> up,low;
    for (int i=0; i<t_dim; ++i) {
      low[i]=double(is[i])/double(ns[i]);
      up[i]=double(is[i]+1)/double(ns[i]);
    }
    Parni_Channel<t_dim> nch;
    nch.SetLower(low);
    nch.SetUpper(up);
    std::cout<<" add: "<<nch;
    m_channels.push_back(nch);

    // increment
    for (int i=0; i<t_dim; ++i) {
      ++is[i];
      if (is[i]<ns[i]) break;
      is[i]=0;
      if (i==t_dim-1) complete=true;
    }
  }
}

template <int t_dim>
void Parni<t_dim>::SplittChannel(const Parni_Channel<t_dim> & mother,int n)
{
  if (n==t_dim) {
    m_channels.push_back(mother);
    std::cout<<"Created 1 new channel"<<std::endl;
    return;
  }  

  double middle=0.5*(mother.Lower()[n]+mother.Upper()[n]);
  double alpha=mother.Alpha()/2.;

  Parni_Channel<t_dim> left,right;
  N_Tuple<t_dim> left_cut,right_cut;
  left_cut=mother.Upper();
  left_cut[n]=middle;
  right_cut=mother.Lower();
  right_cut[n]=middle;
  left.SetLower(mother.Lower());
  left.SetUpper(left_cut);
  left.SetAlpha(alpha);
  right.SetLower(right_cut);
  right.SetUpper(mother.Upper());
  right.SetAlpha(alpha);

  SplittChannel(left,n+1);
  SplittChannel(right,n+1);
}

template <int t_dim>
void Parni<t_dim>::SplittChannel(Channel_Iterator it)
{
  N_Tuple<t_dim> low,middle,up;
  low=it->Lower();
  up =it->Upper();
  for (int i=0;i<t_dim;++i) middle[i]=0.5*(low[i]+up[i]);

  double alpha=it->Alpha()/double(t_dim);
  // add two new channels per dimension


  if ((m_mode&1)==1) {
    m_channels.erase(it);
    std::cout<<"deleted 1 old channels"<<std::endl;


    // t_dim division
    for (int i=0; i<t_dim;++i) {
      Parni_Channel<t_dim> left,right;
      N_Tuple<t_dim> left_cut,right_cut;
      left_cut=up;
      left_cut[i]=middle[i];
      right_cut=low;
      right_cut[i]=middle[i];
      left.SetLower(low);
      left.SetUpper(left_cut);
      left.SetAlpha(alpha);
      right.SetLower(right_cut);
      right.SetUpper(up);
      right.SetAlpha(alpha);
      m_channels.push_back(left);
      m_channels.push_back(right);
      std::cout<<"Created 2 new channels"<<std::endl;
    }
  }
  else {
    Parni_Channel<t_dim> mother=(*it);
    mother.SetAlpha(it->Alpha()*2.);

    m_channels.erase(it);
    std::cout<<"deleted 1 old channels"<<std::endl;
    // 2^t_dim
    SplittChannel(mother,0);
  }
}

template <int t_dim>
void Parni<t_dim>::Merge(int nmin) 
{
  if ((m_mode&1)==0) return;

  std::cout<<" would merge now "<<std::endl;
  /*
  // check for equal channels
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    for (Channel_Iterator it2=it;it2!=m_channels.end();++it2) {
      if (it==it2) continue;
      if (*it == *it2) {
	std::cout<<" Merging :"<<std::endl;
	std::cout<<(*it);
	std::cout<<(*it2);

	it->SetAlpha(it->Alpha()+it2->Alpha());
	it2=m_channels.erase(it2);
	it=m_channels.begin();
      }
    }
  }

  if (m_channels.size()>nmin) {
    // check for unimportant channels
    std::cout<<" Merging unimportant channels "<<std::endl;
    double min_alpha=m_channels.front().Alpha();
    Channel_Iterator wina=m_channels.begin();
    int id=0;
    for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
      it->SetId(id++);
      if (it->Alpha()<min_alpha) {
	std::cout<<" found a "<<it->Alpha()<<std::endl;
	min_alpha=it->Alpha();
	wina=it;
      }
    }
    N_Tuple<t_dim> low=wina->Lower();
    N_Tuple<t_dim> up =wina->Upper();
    double alpha=wina->Alpha();
    std::cout<<" Merging :"<<std::endl;
    std::cout<<(*wina);
    m_channels.erase(wina);

    Channel_Iterator winb=m_channels.begin();
    min_alpha=m_channels.front().Alpha();
    for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
      if (it->Alpha()<min_alpha) {
	std::cout<<" found b "<<it->Alpha()<<std::endl;
	min_alpha=it->Alpha();
	winb=it;
      }
    }

    {
      std::cout<<(*winb);
    }
    alpha+=winb->Alpha();
    for (int i=0;i<t_dim;++i) if (winb->Lower()[i]<low[i]) low[i]=winb->Lower()[i];
    for (int i=0;i<t_dim;++i) if (winb->Upper()[i]>up[i])  up[i]=winb->Upper()[i];
    
    m_channels.erase(winb);

    Parni_Channel<t_dim> nch;
    nch.SetLower(low);
    nch.SetUpper(up);
    nch.SetAlpha(alpha);
    std::cout<<nch;
    (*wina)=nch;

    m_channels.push_back(nch);

  }
  */
}

template <int t_dim>
void Parni<t_dim>::WriteOut(const std::string & pid)
{
  std::ofstream ofile(pid.c_str());

  ofile<<m_channels.size()<<" "<<m_name<<std::endl;
  ofile.precision(12);
  for (Channel_Iterator it=m_channels.begin();it!=m_channels.end();++it) {
    ofile<<(*it);
  }
  ofile.close();
}

template <int t_dim>
void Parni<t_dim>::ReadIn(const std::string & pid)
{
  m_channels.clear();

  std::ifstream ifile(pid.c_str());
  std::string buffer;
  int number;
  ifile>>number;
  getline(ifile,buffer);
  std::cout<<buffer<<" "<<number<<std::endl;
  for (int i=0;i<number;++i) {
    getline(ifile,buffer);
    size_t  a=4;
    size_t  b=buffer.find("(");
    size_t  c;
    std::string salpha=buffer.substr(a,b-a);
    buffer=buffer.substr(b);
    char * err;
    double alpha=strtod(salpha.c_str(),&err);
    N_Tuple<t_dim> low,up;
    for (int i=0;i<t_dim;++i) {
      a=buffer.find("(");
      b=buffer.find(",");
      c=buffer.find(")");
      std::cout<<"#"<<buffer.substr(a+1,b-a-1)<<"#"<<buffer.substr(b+1,c-b-1)<<"#"<<std::endl;
      low[i]=strtod(buffer.substr(a+1,b-a-1).c_str(),&err);
      up[i]=strtod(buffer.substr(b+1,c-b-1).c_str(),&err);
      buffer=buffer.substr(c+1);
    }
    Parni_Channel<t_dim> nch;
    nch.SetLower(low);
    nch.SetUpper(up);
    nch.SetAlpha(alpha);
    std::cout<<nch;
    m_channels.push_back(nch);
  }
  std::cout<<" Parni<"<<t_dim<<">::ReadIn finished with "<<m_channels.size()<<" channels "<<std::endl;
}


// ----------------------------------------



template class N_Tuple<1>;
template class N_Tuple<2>;

template class Parni<1>;
template class Parni<2>;
