#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/PreSample_Multi_Channel.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

double fakt(const double &n) 
{
  if (n<=0.0) return 1.0;
  return n*fakt(n-1.0);
}

class Order_OType {
public:
  int operator()(const Multi_Channel* a, const Multi_Channel* b) {
    //    if "a < b" return 1  else 0;
    if (a->OType()<b->OType()) return 1;
    return 0;
  }
};

PreSample_Multi_Channel::PreSample_Multi_Channel
(Phase_Space_Handler *const psh,Color_Integrator *const colint):
  Multi_Channel("Color Sample VHAAG"), p_colint(colint)
{ 
  m_subchannels.clear();
  m_alpha.clear();
  m_palpha.clear();
  p_ran=new double[1];
}

PreSample_Multi_Channel::~PreSample_Multi_Channel()
{
  DropAllChannels();
  if (s1) { delete[] s1; s1 = 0; }
  if (s2) { delete[] s2; s2 = 0; }
  delete[] p_ran; 
}

void PreSample_Multi_Channel::AddMC(Multi_Channel *mc)
{ 
  m_subchannels.push_back(mc);
  double w=1./m_subchannels.size();
  for (size_t i=0;i<m_alpha.size();i++) m_alpha[i]*=1.-w;
  for (size_t i=0;i<m_alpha.size();i++) m_palpha[i]*=1.-w;
  m_alpha.push_back(w);
  m_palpha.push_back(w);
  m_calpha = w;

  nin = (*mc)[0]->Nin();
  nout = (*mc)[0]->Nout();
  w    = fakt(nout)*(16*pow(3.,nout-2*(nin+nout)));
  if((1<<((nout+1)/2))==(*mc)[0]->OType()) w/=2.;
  m_multi.push_back(w);
}

std::vector<Single_Channel*> PreSample_Multi_Channel::ExtractChannels()
{
  std::vector<Single_Channel*> chs;
  for(size_t i(0);i<m_subchannels.size();++i) {
    for(size_t j(0);j<m_subchannels[i]->Number();++j)
      m_subchannels[i]->Channel(j)->SetAlpha(m_alpha[i]);
    msg_Tracking()<<"extract weight "<<m_subchannels[i]->Name()
		  <<" "<<m_alpha[i]<<"\n";
    chs.insert(chs.end(),m_subchannels[i]->Channels().begin(),
	       m_subchannels[i]->Channels().end());
    m_subchannels[i]->Channels().clear();
  }
  return chs;
}

void PreSample_Multi_Channel::DropAllChannels(const bool del)
{
  for(size_t i=m_subchannels.size();i>0;i--) {
    if (m_subchannels[i-1]) delete m_subchannels[i-1];
  }
  m_subchannels.clear();
}

void PreSample_Multi_Channel::Reset() 
{
  if (s1==0) s1 =  new double[m_subchannels.size()];
  if (s2==0) s2 =  new double[m_subchannels.size()];
  if (!m_readin) {
    s1xmin     = 1.e32;
    n_points   = 0;  
    n_contrib  = 0;
  }

  msg_Tracking()<<"Channels before sort: "<<m_subchannels.size()<<endl;
   std::stable_sort(m_subchannels.begin(),m_subchannels.end(),Order_OType());
  msg_Tracking()<<"Channels after sort: "<<m_subchannels.size()<<endl;

  for (size_t i=0;i<m_subchannels.size();i++) m_subchannels[i]->Reset();

  msg_Tracking()<<"Channels for "<<name<<endl
		<<"----------------- "<<n_points<<" --------------------"<<endl;
  for(size_t i=0;i<m_subchannels.size();i++) {
    msg_Tracking()<<m_subchannels[i]->OType()<<" "<<(*m_subchannels[i])[0]->Name()<<" "<<i<<" : "<<m_subchannels[i]->Name()<<"  : "<<m_alpha[i]<<endl;
  }
  msg_Tracking()<<"----------------- "<<n_points<<" --------------------"<<endl;
  m_readin=false;

  int vs=60*(int)pow(2.,(*m_subchannels[0])[0]->Nout()-2);
  p_vegas=new Vegas(1,vs,std::string("SMC"));
  //    p_vegas->SetAutoOptimize(500);
}

void PreSample_Multi_Channel::ResetOpt() 
{
  n_points = 0;
  for (size_t i=0;i<m_subchannels.size();i++) m_subchannels[i]->ResetOpt();
}        

void PreSample_Multi_Channel::Optimize(double error)
{
  msg_Tracking()<<"Optimize PreSample_Multi_Channel : "<<name<<endl; 

//    for (size_t i=0;i<m_subchannels.size();i++) 
//      if (m_subchannels[i]->ValidN()>m_subchannels[i]->Number()*100)
//        m_subchannels[i]->Optimize(error);
  if (m_fixalpha) return;

  double aptot = 0.,s1sum=0.,ce=0.84;
  int tp=-1;
  double s(0.0), s2(0.0), n(0.0);
  map<int,double> smap;
  for (size_t i=0;i<m_subchannels.size();i++) {
    if (m_subchannels[i]->OType()!=tp) {
      if (tp>0)  {
	double s1=s2/n-sqr(s/n);
	msg_Tracking()<<"Type "<<tp<<" :  (var/<x>/<x^2>/n) "<<s1
		      <<" / "<<s/n<<" / "<<s2/n<<" / "<<n<<endl;
	smap[tp]=s1;
      }
      tp=m_subchannels[i]->OType();
      s=0.;s2=0.;n=0.;
    }
    s+= m_subchannels[i]->Sum();
    s2+=m_subchannels[i]->Sum2();
    n+= m_subchannels[i]->NN();
  }
  double s1l=s2/n-sqr(s/n);
  smap[tp]=s1l;
  msg_Tracking()<<"Type "<<tp<<" :  (var/<x>/<x^2>/n) "<<s1l
		<<" / "<<s/n<<" / "<<s2/n<<" / "<<n<<endl;

  for (size_t i=0;i<m_subchannels.size();i++) {
    s1[i] = smap[m_subchannels[i]->OType()];
    aptot += pow(m_alpha[i],ce)*pow(sqrt(s1[i]),1.-ce);
    s1sum+=sqrt(s1[i]);
  }
  
  double s1x = 0.;  
  for (size_t i=0;i<m_subchannels.size();i++) {
    if (s1[i]>0.) {
      if (dabs(aptot-sqrt(s1[i]))>s1x) s1x = dabs(aptot-sqrt(s1[i]));
      m_alpha[i]=pow(m_alpha[i],ce)*pow(sqrt(s1[i]),1.-ce)/aptot;
//       m_alpha[i]=sqrt(m_alpha[i]*sqrt(s1[i]))/aptot;
    }
  }

  aptot = 0.;
  for (size_t i=0;i<m_subchannels.size();i++) 
    aptot+=1./sqrt(Max(0.5,(double)m_subchannels[i]->ValidN()));
  for (size_t i=0;i<m_subchannels.size();i++) 
    m_palpha[i]=1./sqrt(Max(0.5,(double)m_subchannels[i]->ValidN()))/aptot;

//   if (s1x<s1xmin) {
//     s1xmin = s1x;
//     for (i=0;i<channels.size();i++) channels[i]->SetAlphaSave(channels[i]->Alpha());
//   }  

  msg_Tracking()<<"New weights for : "<<name<<endl
		<<"----------------- "<<n_points<<" ----------------"<<endl;
  for (size_t i=0;i<m_subchannels.size();i++) {
    msg_Tracking()<<i<<" channel "<<m_subchannels[i]->Name()<<", "<<m_subchannels[i]->NN()
		  <<", "<<m_subchannels[i]->N()<<", "<<m_subchannels[i]->ValidN()<<" : "
		  <<m_alpha[i]<<"  ("<<sqrt(s1[i])/s1sum<<")"<<endl;
  }
  msg_Tracking()<<"n,n_contrib : "<<n_points<<", "<<n_contrib<<endl
		<<"-----------------------------------------------"<<endl;
  m_optcnt++;
   for (size_t i=0;i<m_subchannels.size();i++)m_subchannels[i]->ResetCnt();
}

void PreSample_Multi_Channel::EndOptimize(double error)
{ 
  for (size_t i=0;i<m_subchannels.size();i++) 
    m_subchannels[i]->EndOptimize(error);

  msg_Tracking()<<"Best weights:-------------------------------"<<endl;
  for (size_t i=0;i<channels.size();i++) {
    msg_Tracking()<<i<<" channel "<<m_subchannels[i]->Name()<<", "<<m_subchannels[i]->N()
		  <<" : "<<m_alpha[i]<<endl;
  }
  msg_Tracking()<<"n,n_contrib : "<<n_points<<", "<<n_contrib<<endl
		<<"-------------------------------------------"<<endl;
}

bool PreSample_Multi_Channel::OptimizationFinished()
{
  for (size_t i=0;i<m_subchannels.size();i++) 
    if (!m_subchannels[i]->OptimizationFinished()) return false;
  return true;
}

void PreSample_Multi_Channel::AddPoint(double value)
{ 
  if (value!=0.) n_contrib++;
  n_points++;
  m_subchannels[m_lastdice]->AddPoint(value);
}

void PreSample_Multi_Channel::GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts)
{ 
  if(m_subchannels.size()==1) {
    m_subchannels[0]->GeneratePoint(p,cuts);
    m_lastdice = 0;
    return;
  }  

  if(0) {  
    m_fixalpha=true;
    m_lastdice=351;
    m_alpha[m_lastdice]=1.; 
    m_subchannels[m_lastdice]->GeneratePoint(p,cuts);
    return;
  }  

  double dice = ran.Get();
  double sum = 0;
  for (size_t i=0;;++i) {
    if (i==m_subchannels.size()) {
      dice = ran.Get();
      i   = 0;
      sum = 0.;
    }
    sum += m_palpha[i];//m_calpha;
    if (sum>dice) {
      m_subchannels[i]->GeneratePoint(p,cuts);
      m_lastdice = i;
      break;
    }
  }  
  msg_Debugging()<<"selected "<<m_subchannels[m_lastdice]->Name()<<"\n";
  p_colint->GenerateType(m_lastdice,true);  
}

void PreSample_Multi_Channel::GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  m_subchannels[m_lastdice]->GenerateWeight(p,cuts);
  m_weight=m_subchannels[m_lastdice]->Weight()/(m_palpha[m_lastdice]/m_multi[m_lastdice]);
//   m_weight=m_subchannels[m_lastdice]->Weight()/(m_calpha/m_multi[m_lastdice]);
}


void PreSample_Multi_Channel::Print() {
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<"----------------------------------------------"<<endl
		      <<"PreSample_Multi_Channel with "<<m_subchannels.size()<<" channels."<<endl;
  for (size_t i=0;i<m_subchannels.size();i++) 
    msg_Out()<<"  "<<m_subchannels[i]->Name()<<" : "<<m_alpha[i]<<endl;
  msg_Out()<<"++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  for (size_t i=0;i<m_subchannels.size();i++) m_subchannels[i]->Print();
  msg_Out()<<"----------------------------------------------"<<endl;
}                 

void PreSample_Multi_Channel::WriteOut(std::string pid)
{
  ofstream ofile;
  ofile.open(pid.c_str());
  ofile.precision(12);
  ofile<<m_subchannels.size()<<" "<<name<<" "<<n_points<<" "<<n_contrib<<" "
       <<s1xmin<<" "<<m_optcnt<<endl;
  for (size_t i=0;i<m_subchannels.size();i++) 
    ofile<<m_alpha[i]<<" "<<m_palpha[i]<<std::endl;
  ofile.close();
  for (size_t i=0;i<m_subchannels.size();i++) {
    string spid=pid+string("_")+m_subchannels[i]->Name();
    m_subchannels[i]->WriteOut(spid);
  }                 
}                 

bool PreSample_Multi_Channel::ReadIn(std::string pID) {
  ifstream ifile;
  ifile.open(pID.c_str());
  if (ifile.bad()) return false;
  size_t      size;
  std::string name;
  ifile>>size>>name;
  if (( size != m_subchannels.size()) || ( name != name) ) {
    msg_Error()<<"Error in PreSample_Multi_Channel::ReadIn("<<pID<<")"<<endl 
	       <<"  PreSample_Multi_Channel file did not coincide with actual PreSample_Multi_Channel: "<<endl
	       <<"  "<<size<<" vs. "<<channels.size()<<" and "
	       <<"  "<<name<<" vs. "<<name<<endl;
    return 0;
  }
  m_readin=true;
  ifile>>n_points>>n_contrib>>s1xmin>>m_optcnt;

  for (size_t i=0;i<m_subchannels.size();i++) {
    ifile>>m_alpha[i]>>m_palpha[i];
  }
  ifile.close();
  for (size_t i=0;i<channels.size();i++) {
    string spid=pID+string("_")+m_subchannels[i]->Name();
    m_subchannels[i]->ReadIn(spid);
  }
  return 1;
}

int PreSample_Multi_Channel::HandicapFactor() 
{ return 10; }
