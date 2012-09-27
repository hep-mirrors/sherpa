#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Poincare.H"
#include <algorithm>
#ifdef USING__MPI
#include "mpi.h"
#endif

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Multi_Channel::Multi_Channel(string _name,int id) : 
  fl(NULL), m_id(id), s1(NULL), s2(NULL), m_readin(false), m_fixalpha(false),
  m_minalpha(0.0)
{
  string help;
  int    pos;
  for (;;) {
    pos  = _name.find(" ");
    if (pos==-1) break;
    help   = _name;
    _name  = help.substr(0,pos) + help.substr(pos+1); 
  }
  name     = _name;
  n_points = n_contrib = 0;
  mn_points = mn_contrib = 0;
  m_lastdice = -1;
  m_optcnt = 0;
  m_pol = 250.;
  m_otype = 0;
}

Multi_Channel::~Multi_Channel() 
{
  DropAllChannels();
  if (s1) { delete[] s1; s1 = 0; }
  if (s2) { delete[] s2; s2 = 0; }
}

void Multi_Channel::Add(Single_Channel * Ch) { 
  channels.push_back(Ch);
  m_otype = m_otype|Ch->OType();
}

size_t Multi_Channel::NChannels() const
{
  size_t nch(0);
  for (size_t i(0);i<channels.size();++i) nch+=channels[i]->NChannels();
  return nch;
} 

Single_Channel * Multi_Channel::Channel(int i) { 
  if ((i<0) || (i>=(int)channels.size())) {
    msg_Error()<<"Multi_Channel::Channel("<<i<<") out of bounds :"
	       <<" 0 < "<<i<<" < "<<channels.size()<<endl;
    return 0;
  }
  return channels[i]; 
}

void Multi_Channel::DropChannel(int i) 
{
  if ((i<0) || (i>(int)channels.size())) {
    msg_Error()<<"Multi_Channel::DropChannel("<<i<<") out of bounds :"
	       <<" 0 < "<<i<<" < "<<channels.size()<<endl;
    return;
  }
  if (channels[i]) delete channels[i];
  for (size_t j=i;j<channels.size()-1;j++) channels[j] = channels[j+1];
  channels.pop_back();
}

void Multi_Channel::DropAllChannels(const bool del)
{
  while (channels.size()) {
    if (del) delete channels.back();
    channels.pop_back();
  }
}

void Multi_Channel::Reset() 
{
  if (channels.size()==0) {
    if (s1!=NULL) delete[] s1; s1=NULL;
    if (s2!=NULL) delete[] s2; s2=NULL;
    return;
  }
  if (s1!=NULL) delete[] s1;
  s1 =  new double[channels.size()];
  if (s2!=NULL) delete[] s2;
  s2 =  new double[channels.size()];
  if (!m_readin) {
    s1xmin     = 1.e32;
    n_points   = 0;  
    n_contrib  = 0;
    mn_points = mn_contrib = 0;
  }
  msg_Tracking()<<"Channels for "<<name<<endl
		<<"----------------- "<<n_points<<" --------------------"<<endl;
  for(size_t i=0;i<channels.size();i++) {
    if (!m_readin) channels[i]->Reset(1./channels.size());
    msg_Tracking()<<" "<<i<<" : "<<channels[i]->Name()<<"  : "<<channels[i]->Alpha()<<endl;
  }
  msg_Tracking()<<"----------------- "<<n_points<<" --------------------"<<endl;
  m_readin=false;
}

void Multi_Channel::ResetOpt() 
{
  n_points = 0;
}        

void Multi_Channel::ResetCnt() 
{
}        

class Order_Weight {
public:
  bool operator()(Single_Channel* c1,Single_Channel* c2)
  { 
    return c1->Alpha()>c2->Alpha();
  }
};

void Multi_Channel::MPISync()
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    int rank=MPI::COMM_WORLD.Get_rank();
    int cn=3*channels.size()+2;
    double *values = new double[cn];
    if (rank==0) {
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	MPI::COMM_WORLD.Recv(values,cn,MPI::DOUBLE,MPI::ANY_SOURCE,tag);
	for (size_t i=0;i<channels.size();++i) {
	  channels[i]->AddMPIVars(values[i],
				  values[channels.size()+i],
				  values[2*channels.size()+i]);
	}
	mn_points+=values[cn-2];
	mn_contrib+=values[cn-1];
      }
      for (size_t i=0;i<channels.size();++i) {
	values[i]=channels[i]->MRes1();
	values[channels.size()+i]=channels[i]->MRes2();
	values[2*channels.size()+i]=channels[i]->MRes3();
      }
      values[cn-2]=mn_points;
      values[cn-1]=mn_contrib;
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	MPI::COMM_WORLD.Send(values,cn,MPI::DOUBLE,tag,size+tag);
      }
    }
    else {
      for (size_t i=0;i<channels.size();++i) {
	values[i]=channels[i]->MRes1();
	values[channels.size()+i]=channels[i]->MRes2();
	values[2*channels.size()+i]=channels[i]->MRes3();
      }
      values[cn-2]=mn_points;
      values[cn-1]=mn_contrib;
      MPI::COMM_WORLD.Send(values,cn,MPI::DOUBLE,0,rank);
      MPI::COMM_WORLD.Recv(values,cn,MPI::DOUBLE,0,size+rank);
      for (size_t i=0;i<channels.size();++i) {
	channels[i]->SetMPIVars(values[i],
				values[channels.size()+i],
				values[2*channels.size()+i]);
      }
      mn_points=values[cn-2];
      mn_contrib=values[cn-1];
    }
    delete [] values;
  }
  for (size_t i=0;i<channels.size();++i) {
    channels[i]->CopyMPIValues();
    channels[i]->MPISync();
  }
  n_points+=mn_points;
  n_contrib+=mn_contrib;
  mn_points=mn_contrib=0.0;
#endif
}

void Multi_Channel::Optimize(double error)
{
  if (m_fixalpha) return;
  msg_Tracking()<<"Optimize Multi_Channel : "<<name<<endl; 

  double aptot = 0.;
  size_t i;
  for (i=0;i<channels.size();i++) {
    s1[i]  = channels[i]->Res1()/n_points;
    s2[i]  = sqrt(channels[i]->Res2()-
		  channels[i]->Res3()/(n_points-1. ))/n_points;
    aptot += channels[i]->Alpha()*sqrt(s1[i]);
  }
  
  double s1x = 0.;  
  for (i=0;i<channels.size();i++) {
    if (channels[i]->Alpha()>0.) {
      if (dabs(aptot-sqrt(s1[i]))>s1x) s1x = dabs(aptot-sqrt(s1[i]));
      if (channels.size()>1) {
      channels[i]->SetAlpha(channels[i]->Alpha() * sqrt(s1[i])/aptot);
      if (channels[i]->Alpha() < Min(1.e-4,1.e-3/(double)channels.size()) ) channels[i]->SetAlpha(m_minalpha);
      }
    }
  }
  double norm = 0;
  for (i=0;i<channels.size();i++) norm += channels[i]->Alpha();
  for (i=0;i<channels.size();i++) channels[i]->SetAlpha(channels[i]->Alpha() / norm);

    if((m_optcnt>4 || channels.size()==1)&& m_optcnt<20)    
      for (i=0;i<channels.size();i++) if (channels[i]->Alpha()>0.01) channels[i]->Optimize();
    if (m_optcnt==20){
      for (i=0;i<channels.size();i++) if (channels[i]->Alpha()>0.) channels[i]->EndOptimize();
      s1xmin     = 1.e32;
    }

  if (s1x<s1xmin) {
    s1xmin = s1x;
    for (i=0;i<channels.size();i++) channels[i]->SetAlphaSave(channels[i]->Alpha());
  }  

  for(i=0;i<channels.size();i++) channels[i]->ResetOpt();
  msg_Tracking()<<"New weights for : "<<name<<endl
		<<"----------------- "<<n_points<<" ----------------"<<endl;
  for (i=0;i<channels.size();i++) {
    if (channels[i]->Alpha() > 0) {
      msg_Tracking()<<i<<" channel "<<channels[i]->Name()<<", "<<channels[i]->N()<<" : "
		    <<channels[i]->Alpha()<<" -> "<<channels[i]->AlphaSave()<<endl;
    }
  }
  msg_Tracking()<<"S1X: "<<s1x<<" -> "<<s1xmin<<endl
 		<<"n,n_contrib : "<<n_points<<", "<<n_contrib<<endl
		<<"-----------------------------------------------"<<endl;
  m_optcnt++;
  m_best=channels;
  std::sort(m_best.begin(),m_best.end(),Order_Weight());
  m_best.resize(2);
}

void Multi_Channel::EndOptimize(double error)
{
  size_t i;

  for (i=0;i<channels.size();i++) {
    channels[i]->SetAlpha(channels[i]->AlphaSave());
    if (channels[i]->Alpha() < Min(1.e-4,1.e-2/(double)channels.size())) channels[i]->SetAlpha(m_minalpha);
  }
  double norm = 0;
  for (i=0;i<channels.size();i++) norm += channels[i]->Alpha();
  for (i=0;i<channels.size();i++) channels[i]->SetAlpha(channels[i]->Alpha() / norm);
  for (i=0;i<channels.size();i++) if (channels[i]->Alpha()>0.) channels[i]->EndOptimize();

  msg_Tracking()<<"Best weights:-------------------------------"<<endl;
  for (i=0;i<channels.size();i++) {
    if (channels[i]->Alpha() > 0) {
      msg_Tracking()<<i<<" channel "<<channels[i]->Name()<<", "<<channels[i]->N()
		    <<" : "<<channels[i]->Alpha()<<endl;
    }
  }
  msg_Tracking()<<"S1X: "<<s1xmin<<endl
 		<<"n,n_contrib : "<<n_points<<", "<<n_contrib<<endl
		<<"-------------------------------------------"<<endl;
}

bool Multi_Channel::OptimizationFinished()
{
  for (size_t i=0;i<channels.size();i++) if (!channels[i]->OptimizationFinished()) return false;
  return true;
}


void Multi_Channel::AddPoint(double value)
{
  //if (!ATOOLS::IsZero(value)) n_contrib++;
#ifdef USING__MPI
  if (value!=0.) mn_contrib++;
#else
  if (value!=0.) n_contrib++;
#endif
  //   if (value!=0.) PRINT_INFO(Name()<<" "<<value<<" "<<n_contrib);
#ifdef USING__MPI
  mn_points++;
#else
  n_points++;
#endif
  double var;
  for (size_t i=0;i<channels.size();i++) {
    if (value!=0.) {
      if (channels[i]->Weight()!=0) 
	var = sqr(value)*m_weight/channels[i]->Weight();
      else var = 0.;
#ifdef USING__MPI
      channels[i]->AddMPIVars(var,sqr(var),
			      sqr(channels[i]->MRes1())
			      -channels[i]->MRes2()
			      -channels[i]->MRes3());
#else
      channels[i]->SetRes1(channels[i]->Res1() + var);
      channels[i]->SetRes2(channels[i]->Res2() + sqr(var));
      channels[i]->SetRes3(sqr(channels[i]->Res1())-channels[i]->Res2());
#endif
    }
  }
  if (m_lastdice>=0) Channel(m_lastdice)->AddPoint(value);
}



void Multi_Channel::GenerateWeight(Vec4D * p,Cut_Data * cuts)
{
  Vec4D_Vector pp(p,&p[nin+nout]);
  if (nin==2) {
    Poincare cms(pp[0]+pp[1]);
    for (int i(0);i<nin+nout;++i) cms.Boost(pp[i]);
  }
  if (channels.size()==1) {
    channels[0]->GenerateWeight(&pp.front(),cuts);
    if (channels[0]->Weight()!=0) m_weight = channels[0]->Weight();
    return;
  }
  m_weight = 0.;
  for (size_t i=0; i<channels.size(); ++i) {
    if (channels[i]->Alpha() > 0.) {
      channels[i]->GenerateWeight(&pp.front(),cuts);
      if (!(channels[i]->Weight()>0) && 
	  !(channels[i]->Weight()<0) && (channels[i]->Weight()!=0)) {
	msg_Error()<<"Multi_Channel::GenerateWeight(..): ("<<this->name
		   <<"): Channel "<<i<<" ("<<channels[i]<<") produces a nan!"<<endl;
      }
      if (channels[i]->Weight()!=0) 
	m_weight += channels[i]->Alpha()/channels[i]->Weight();
    }
  }
  if (m_weight!=0) m_weight = 1./m_weight;
}


void Multi_Channel::GeneratePoint(Vec4D *p,Cut_Data * cuts)
{
  Poincare cms(p[0]+p[1]);
  if (nin==2) for (int i(0);i<nin;++i) cms.Boost(p[i]);
  for(size_t i=0;i<channels.size();i++) channels[i]->SetWeight(0.);
  if(channels.size()==1) {
    channels[0]->GeneratePoint(p,cuts);
    if (nin==2) for (int i(0);i<nin+nout;++i) cms.BoostBack(p[i]);
    m_lastdice = 0;
    return;
  }  
  double rn  = ran->Get();
  double sum = 0;
  for (size_t i=0;;++i) {
    if (i==channels.size()) {
      rn  = ran->Get();
      i   = 0;
      sum = 0.;
    }
    sum += channels[i]->Alpha();
    if (sum>rn) {
      channels[i]->GeneratePoint(p,cuts);
      if (nin==2) for (int i(0);i<nin+nout;++i) cms.BoostBack(p[i]);
      m_lastdice = i;
      break;
    }
  }  
}

void Multi_Channel::GeneratePoint(Info_Key &spkey,Info_Key &ykey,int mode) 
{
  for(size_t i=0;i<channels.size();++i) channels[i]->SetWeight(0.);
  double disc=ran->Get();
  double sum=0.;
  for (size_t n=0;n<channels.size();++n) {
    sum+=channels[n]->Alpha();
    if (sum>disc) {
      for (size_t i=0;i<2;++i) rans[i]=ran->Get();
      channels[n]->GeneratePoint(spkey,ykey,rans,mode);
      m_lastdice = n;
      return;
    }
  }  
  msg_Error()<<"Multi_Channel::GeneratePoint(..): IS case ("<<this
	     <<") No channel selected. \n"
	     <<"   disc = "<<disc<<", sum = "<<sum<<std::endl;
  abort();
}

void Multi_Channel::GenerateWeight(int mode=0)
{
  if (channels.size()==1) {
    channels[0]->GenerateWeight(mode);
    if (channels[0]->Weight()!=0) m_weight = channels[0]->Weight();
    return;
  }
  m_weight = 0.;
  for (size_t i=0;i<channels.size();++i) {
    if (channels[i]->Alpha()>0.) {
      channels[i]->GenerateWeight(mode);
      if (!(channels[i]->Weight()>0)&&
	  !(channels[i]->Weight()<0)&&(channels[i]->Weight()!=0)) {
	msg_Error()<<"Multi_Channel::GenerateWeight(): ("<<this->name
		   <<"): Channel "<<i<<" ("<<channels[i]<<") produces a nan!"<<endl;
      }
      if (channels[i]->Weight()!=0) 
	m_weight += channels[i]->Alpha()/channels[i]->Weight();
    }
  }
  if (m_weight!=0) m_weight=1./m_weight;
}

void Multi_Channel::ISRInfo(int i,int & type,double & mass,double & width) 
{
  channels[i]->ISRInfo(type,mass,width);
  return;
}

void Multi_Channel::ISRInfo
(std::vector<int> &ts,std::vector<double> &ms,std::vector<double> &ws) const
{
  for (size_t i=0;i<channels.size();++i) channels[i]->ISRInfo(ts,ms,ws);
}

void Multi_Channel::Print() {
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<"----------------------------------------------"<<endl
		      <<"Multi_Channel with "<<channels.size()<<" channels."<<endl;
  for (size_t i=0;i<channels.size();i++) 
    msg_Out()<<"  "<<channels[i]->Name()<<" : "<<channels[i]->Alpha()<<endl;
  msg_Out()<<"----------------------------------------------"<<endl;
}                 


void Multi_Channel::WriteOut(std::string pID) 
{
  ofstream ofile;
  ofile.open(pID.c_str());
  ofile.precision(12);
  ofile<<channels.size()<<" "<<name<<" "<<n_points<<" "<<n_contrib<<" "
       <<s1xmin<<" "<<m_optcnt<<endl;
//        <<m_result<<" "<<m_result2<<" "<<s1xmin<<" "
//        <<m_sresult<<" "<<m_sresult2<<" "<<m_ssigma2<<" "<<n_spoints<<" "<<m_optcnt<<endl;
  for (size_t i=0;i<channels.size();i++) 
    ofile<<channels[i]->Name()<<" "<<channels[i]->N()<<" "
	 <<channels[i]->Alpha()<<" "<<channels[i]->AlphaSave()<<" "
	 <<channels[i]->Weight()<<" "<<channels[i]->Res1()<<" "
	 <<channels[i]->Res2()<<" "<<channels[i]->Res3()<<std::endl;
  ofile.close();
  for (size_t i=0;i<channels.size();i++) channels[i]->WriteOut(pID);
}

bool Multi_Channel::ReadIn(std::string pID) {
  ifstream ifile;
  ifile.open(pID.c_str());
  if (ifile.bad()) return false;
  size_t      size;
  std::string name;
  long int    points;
  double      alpha, alphasave, weight, res1, res2, res3;
  ifile>>size>>name;
  if (( size != channels.size()) || ( name != name) ) {
    msg_Error()<<"Error in Multi_Channel::ReadIn("<<pID<<")"<<endl 
	       <<"  Multi_Channel file did not coincide with actual Multi_Channel: "<<endl
	       <<"  "<<size<<" vs. "<<channels.size()<<" and "
	       <<"  "<<name<<" vs. "<<name<<endl;
    return 0;
  }
  m_readin=true;
  //   ifile>>n_points>>n_contrib>>m_result>>m_result2>>s1xmin>>m_sresult
  // >>m_sresult2>>m_ssigma2>>n_spoints>>m_optcnt;
  ifile>>n_points>>n_contrib>>s1xmin>>m_optcnt;

  double sum=0;
  for (size_t i=0;i<channels.size();i++) {
    ifile>>name>>points>>alpha>>alphasave>>weight>>res1>>res2>>res3;
    sum+= alpha;
    if (name != channels[i]->Name()) {
      msg_Error()<<"ERROR in "<<METHOD<<" for "<<pID<<")"<<endl 
		 <<"  name of Single_Channel not consistent ("<<i<<")"<<endl
		 <<"  "<<name<<" vs. "<<channels[i]->Name()<<endl;
      return 0;
      if (name.substr(0,name.length()-1)!=
	  channels[i]->Name().substr(0,name.length()-1)) {
	msg_Error()<<"   return 0."<<std::endl;
	return 0;
      }
    }
    channels[i]->SetN(points);
    channels[i]->SetAlpha(alpha);
    channels[i]->SetAlphaSave(alphasave);
    channels[i]->SetWeight(weight);
    channels[i]->SetRes1(res1);
    channels[i]->SetRes2(res2);
    channels[i]->SetRes3(res3);
  }
  ifile.close();
  for (size_t i=0;i<channels.size();i++) channels[i]->ReadIn(pID);
  return 1;
}

std::string Multi_Channel::ChID(int n)
{
  return channels[n]->ChID();
}

void Multi_Channel::SetRange(double *_sprimerange,double *_yrange) 
{
  for (size_t i=0;i<channels.size();++i) channels[i]->SetRange(_sprimerange,_yrange);
}

void Multi_Channel::GetRange() 
{
  for (unsigned int i=0;i<channels.size();i++) channels[i]->GetRange();
}

bool Multi_Channel::Initialize()
{ 
  return true;
}
