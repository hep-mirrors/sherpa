#include "Multi_Channel.H"
#include "Random.H"

using namespace PHASIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

//#define _USE_MPI_
#ifdef _USE_MPI_
#include <mpi++.h>
#endif

Multi_Channel::Multi_Channel(string _name) : fl(NULL), s1(NULL), s2(NULL)
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
}

Multi_Channel::~Multi_Channel() 
{
  DropAllChannels();
  if (s1) { delete[] s1; s1 = 0; }
  if (s2) { delete[] s2; s2 = 0; }
  msg.Debugging()<<"Deleted "<<name<<endl;
}

void Multi_Channel::Add(Single_Channel * Ch) { 
  AORGTOOLS::msg.Debugging()<<"Add "<<Ch->Name()<<" to the multi-channel."<<endl;
  channels.push_back(Ch); 
}

Single_Channel * Multi_Channel::Channel(int i) { 
  if ((i<0) || (i>=channels.size())) {
    msg.Error()<<"Multi_Channel::Channel("<<i<<") out of bounds :";
    msg.Error()<<" 0 < "<<i<<" < "<<channels.size()<<endl;
    return 0;
  }
  return channels[i]; 
}

void Multi_Channel::DropChannel(int i) 
{
  if ((i<0) || (i>channels.size())) {
    msg.Error()<<"Multi_Channel::DropChannel("<<i<<") out of bounds :";
    msg.Error()<<" 0 < "<<i<<" < "<<channels.size()<<endl;
    return;
  }
  if (channels[i]) delete channels[i];
  for (short int j=i;j<channels.size()-1;j++) channels[j] = channels[j+1];
  channels.pop_back();
}

void Multi_Channel::DropAllChannels()
{
  for(int i=channels.size();i>0;i--) {
    if (channels[i-1]) delete channels[i-1];
  }
  channels.clear();
  AORGTOOLS::msg.Debugging()<<"Dropped all channels for Multi_Channel : "<<name<<endl;
}

void Multi_Channel::Reset() 
{
  msg.Debugging()<<"Resetting Multi_Channel : "<<this<<endl;
  if (s1==0) s1 =  new double[channels.size()];
  if (s2==0) s2 =  new double[channels.size()];

  s1xmin     = 1.e32;
  n_points   = 0;  
  n_contrib  = 0;
  m_result     = 0.;
  m_result2    = 0.;

  msg.Tracking()<<"Channels for "<<name<<endl;
  msg.Tracking()<<"----------------- "<<n_points<<" --------------------"<<endl;
  for(short int i=0;i<channels.size();i++) {
    channels[i]->Reset(1./channels.size());
    msg.Tracking()<<" "<<i<<" : "<<channels[i]->Name()<<"  : "<<channels[i]->Alpha()<<endl;
  }
  msg.Tracking()<<"----------------- "<<n_points<<" --------------------"<<endl;
}

void Multi_Channel::ResetOpt() 
{
  n_points = 0;
}        

void Multi_Channel::MPIOptimize(double error)
{
#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  
  cout<<"Process "<<rank<<" in MPIOptimize()."<<endl;

  double * messageblock = new double[1+3*channels.size()];
  double * alp = new double[channels.size()];
  //tag 9 for communication
  if (rank==0) {
    int count = 0;
    while (count<size-1) {
      MPI::COMM_WORLD.Recv(messageblock, 1+3*channels.size(), MPI::DOUBLE, MPI::ANY_SOURCE, 9);
      count++;
      cout<<" received Channel-Data from Knot "<<messageblock[0]<<endl;
      for (short int i=0;i<channels.size();i++) {
	channels[i]->SetRes1(channels[i]->Res1() + messageblock[1+channels.size()+i]);
	channels[i]->SetRes2(channels[i]->Res2() + messageblock[1+2*channels.size()+i]);
	channels[i]->SetRes3(sqr(channels[i]->Res1()) - channels[i]->Res2());
	channels[i]->SetN(channels[i]->N() + int(messageblock[1+i]));
      }
    }
    cout<<"Master: "<<channels[0]->N()<<" data accumulated."<<endl; 
    Optimize(error);
    //broadcast alphai
    for (short int i=0;i<channels.size();i++) alp[i] = channels[i]->Alpha();
    
    for (int i=1;i<size;i++) MPI::COMM_WORLD.Isend(alp, channels.size(), MPI::DOUBLE, i, 9);
  }
  else {
    messageblock[0] = rank;
    for (int i=0;i<channels.size();i++) messageblock[1+i]                   = channels[i]->N();
    for (int i=0;i<channels.size();i++) messageblock[1+channels.size()+i]   = channels[i]->Res1();
    for (int i=0;i<channels.size();i++) messageblock[1+2*channels.size()+i] = channels[i]->Res2();
    MPI::COMM_WORLD.Send(messageblock, 1+3*channels.size(), MPI::DOUBLE, 0, 9);   
    //Waiting for new alpha's
    MPI::COMM_WORLD.Recv(alp, channels.size(), MPI::DOUBLE, 0, 9);
    msg.Out()<<"Slave "<<rank<<" received new alpha."<<endl;

    for (short int i=0;i<channels.size();i++) {
      channels[i]->SetAlpha(alp[i]);
      channels[i]->ResetOpt();
      channels[i]->SetWeight(0.);
    }
  }
  delete[] alp;
  delete[] messageblock;
#else
  cout<<"MPIOptimize called in non-MPI session!!!"<<endl;
#endif

}


void Multi_Channel::Optimize(double error)
{
  msg.Tracking()<<"Optimize Multi_Channel : "<<name<<endl; 

  double aptot = 0.;
  short int i;
  for (i=0;i<channels.size();i++) {
    s1[i]  = channels[i]->Res1()/channels[i]->N();
    s2[i]  = sqrt(channels[i]->Res2()-
		  channels[i]->Res3()/(channels[i]->N()-1. ))/channels[i]->N();
    aptot += channels[i]->Alpha()*sqrt(s1[i]);
  }
  
  double s1x = 0.;  
  for (i=0;i<channels.size();i++) {
    if (dabs(aptot-sqrt(s1[i]))>s1x) s1x = dabs(aptot-sqrt(s1[i]));
    channels[i]->SetAlpha(channels[i]->Alpha() * sqrt(s1[i])/aptot);
    //maximum number of events equal 10^8 assumed
    if (channels[i]->Alpha() < 1.e-8 ) channels[i]->SetAlpha(0.);
    //    if (channels[i]->Alpha() < sqr(error)/channels.size()) channels[i]->SetAlpha(0.);
  }
  double norm = 0;
  for (i=0;i<channels.size();i++) norm += channels[i]->Alpha();
  for (i=0;i<channels.size();i++) channels[i]->SetAlpha(channels[i]->Alpha() / norm);

  if (s1x<s1xmin) {
    s1xmin = s1x;
    for (i=0;i<channels.size();i++) channels[i]->SetAlphaSave(channels[i]->Alpha());
  }  
  for(i=0;i<channels.size();i++) channels[i]->ResetOpt();
  msg.Tracking()<<"New weights for : "<<name<<endl;
  msg.Tracking()<<"----------------- "<<n_points<<" ----------------"<<endl;
  for (i=0;i<channels.size();i++) {
    if (channels[i]->Alpha() > 0) {
      msg.Tracking()<<i<<" channel "<<channels[i]->Name()<<", "<<channels[i]->N()<<" : ";
      msg.Tracking()<<channels[i]->Alpha()<<" -> "<<channels[i]->AlphaSave()<<endl;
    }
  }
  msg.Tracking()<<"S1X: "<<s1x<<" -> "<<s1xmin<<endl;
  msg.Tracking()<<"Variance : "<<Variance()<<endl;
  msg.Tracking()<<"result,result2,n,n_contrib : "<<m_result<<", ";
  msg.Tracking()<<m_result2<<", "<<n_points<<", "<<n_contrib<<endl;
  msg.Tracking()<<"-----------------------------------------------"<<endl;
}

void Multi_Channel::EndOptimize(double error)
{
  short int i;

#ifndef _USE_MPI_

  for (i=0;i<channels.size();i++) {
    channels[i]->SetAlpha(channels[i]->AlphaSave());
    if (channels[i]->Alpha() < 1.e-8 ) channels[i]->SetAlpha(0.);
    //    if (channels[i]->Alpha() < error/channels.size()) channels[i]->SetAlpha(0.);
  }
  double norm = 0;
  for (i=0;i<channels.size();i++) norm += channels[i]->Alpha();
  for (i=0;i<channels.size();i++) channels[i]->SetAlpha(channels[i]->Alpha() / norm);

  msg.Tracking()<<"Best weights:-------------------------------"<<endl;
  for (i=0;i<channels.size();i++) {
    if (channels[i]->Alpha() > 0) {
      msg.Tracking()<<i<<" channel "<<channels[i]->Name()<<", "<<channels[i]->N();
      msg.Tracking()<<" : "<<channels[i]->Alpha()<<endl;
    }
  }
  msg.Tracking()<<"S1X: "<<s1xmin<<endl;
  msg.Tracking()<<"Variance : "<<Variance()<<endl;
  msg.Tracking()<<"result,result2,n,n_contrib : "<<m_result<<", ";
  msg.Tracking()<<m_result2<<", "<<n_points<<", "<<n_contrib<<endl;
  msg.Tracking()<<"-------------------------------------------"<<endl;

#else

//bcast them to all
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  
  cout<<"Process "<<rank<<" in End_Optimize()."<<endl;

  double* alp = new double[channels.size()];
  //tag 9 for communication
  if (rank==0) {
    short int i;
    for (i=0;i<channels.size();i++) channels[i]->SetAlpha(channels[i]->AlphaSave());
    double norm = 0;
    for (i=0;i<channels.size();i++) norm += channels[i]->Alpha();
    for (i=0;i<channels.size();i++) channels[i]->SetAlpha(channels[i]->Alpha() / norm);

    msg.Tracking()<<"Best weights:-------------------------------"<<endl;
    for (i=0;i<channels.size();i++)
      if (channels[i]->Alpha() > 0) {
	msg.Tracking()<<i<<" channel "<<channels[i]->Name()<<" :"<<channels[i]->Alpha()<<endl;
	msg.Tracking()<<"S1X: "<<s1xmin<<endl;
	msg.Tracking()<<"-------------------------------------------"<<endl;
      }
    //broadcast alphai
    for (short int i=0;i<channels.size();i++) alp[i] = channels[i]->Alpha();    
    for (int i=1;i<size;i++) MPI::COMM_WORLD.Isend(alp, channels.size(), MPI::DOUBLE, i, 9);
  }
  else {
    //Waiting for new alpha's
    MPI::COMM_WORLD.Recv(alp, channels.size(), MPI::DOUBLE, 0, 9);
    cout<<"Slave "<<rank<<" received new alpha."<<endl;
    for (short int i=0;i<channels.size();i++) channels[i]->SetAlpha(alp[i]);
  }
  delete[] alp;
#endif
}

void Multi_Channel::AddPoint(double value)
{
  // msg.Debugging()<<"In Multi_Channel::AddPoint("<<value<<")"<<endl;
  if (!AMATOOLS::IsZero(value)) n_contrib++;

  n_points++;
  m_result  += value;
  m_result2 += value*value;
  
  double var;
  for (short int i=0;i<channels.size();i++) {
    if (value!=0.) {
      if (channels[i]->Weight()!=0) 
	var = sqr(value)*m_weight/channels[i]->Weight();
      else var = 0.;

      //Reciprocal weights compared to Berends et al.
      channels[i]->SetRes1(channels[i]->Res1() + var);
      channels[i]->SetRes2(channels[i]->Res2() + sqr(var));
      channels[i]->SetRes3(sqr(channels[i]->Res1())-channels[i]->Res2());
    }
    channels[i]->IncrementN();
  }
}

double Multi_Channel::Variance() {
  double disc = (n_points*m_result2)/((n_points-1)*AMATOOLS::sqr(m_result)) - 1./(n_points-1);
  if (disc>0.) return m_result/n_points * sqrt(disc);
  disc = m_result2/(n_points*(n_points-1)) - AMATOOLS::sqr(m_result/n_points)/(n_points-1);
  if (disc>0.) return sqrt(disc);
  
  AORGTOOLS::msg.Error()<<"Variance yielded a NaN !"<<endl;
  AORGTOOLS::msg.Error()<<"   res,res2 = "<<m_result<<", "<<m_result2;
  AORGTOOLS::msg.Error()<<" after "<<n_points<<" points."<<endl; 
  
  return sqrt(-disc);
};


void Multi_Channel::GenerateWeight(int n,Vec4D* p,Cut_Data * cuts) 
{
  if (channels[n]->Alpha() > 0.) {
    channels[n]->GenerateWeight(p,cuts);
    if (channels[n]->Weight()==0.) m_weight = 0.; 
                              else m_weight = 1./channels[n]->Weight();
  }
  else m_weight = 0.;
}

void Multi_Channel::GenerateWeight(Vec4D * p,Cut_Data * cuts)
{
  m_weight = 0.;
  for (short int i=0; i<channels.size(); ++i) {
    if (channels[i]->Alpha() > 0.) {
      channels[i]->GenerateWeight(p,cuts);
      if (!(channels[i]->Weight()>0) && 
	  !(channels[i]->Weight()<0) && (channels[i]->Weight()!=0)) {
	msg.Error()<<"Channel "<<i<<" produces a nan!"<<endl;
      }
      if (channels[i]->Weight()!=0) 
	m_weight += channels[i]->Alpha()/channels[i]->Weight();
    }
  }
  if (m_weight!=0) m_weight = 1./m_weight;
}


void Multi_Channel::GeneratePoint(int n,Vec4D * p,Cut_Data * cuts,double * ran)
{
  channels[n]->GeneratePoint(p,cuts,ran);
}

void Multi_Channel::GeneratePoint(Vec4D * p,Cut_Data * cuts)
{
  for(short int i=0;i<channels.size();i++) channels[i]->SetWeight(0.);
  if(channels.size()==1) {
    channels[0]->GeneratePoint(p,cuts);
    return;
  }  
  double rn  = ran.Get();
  double sum = 0;
  for (short int i=0;;++i) {
    if (i==channels.size()) {
      rn  = ran.Get();
      i   = 0;
      sum = 0.;
    }
    sum += channels[i]->Alpha();
    if (sum>rn) {
      //       cout<<"Channel number "<<i<<"  rn="<<rn<<" sum="<<sum<<endl;
      channels[i]->GeneratePoint(p,cuts);
      break;
    }
  }  
}

void Multi_Channel::GenerateWeight(int n,Vec4D* p) 
{
  if (channels[n]->Alpha() > 0.) {
    channels[n]->GenerateWeight(p);
    if (channels[n]->Weight()==0.) m_weight = 0.; 
                              else m_weight = 1./channels[n]->Weight();
  }
  else m_weight = 0.;
}


void Multi_Channel::GenerateWeight(Vec4D * p)
{
  m_weight = 0.;
  for (short int i=0; i<channels.size(); ++i) {
    if (channels[i]->Alpha() > 0.) {
      channels[i]->GenerateWeight(p);
      if (!(channels[i]->Weight()>0) && 
	  !(channels[i]->Weight()<0) && (channels[i]->Weight()!=0)) {
	msg.Error()<<"Channel "<<i<<" produces a nan!"<<endl;
      }
      if (channels[i]->Weight()!=0) 
	m_weight += channels[i]->Alpha()/channels[i]->Weight();
    }
  }
  if (!AMATOOLS::IsZero(m_weight)) m_weight = 1./m_weight;
}


void Multi_Channel::GeneratePoint(int n,Vec4D * p,double * rn)
{
  channels[n]->GeneratePoint(p,rn);
}

void Multi_Channel::GeneratePoint(Vec4D * p)
{
  for(short int i=0;i<channels.size();i++) channels[i]->SetWeight(0.);
  double rn  = ran.Get();
  double sum = 0;
  for (short int i=0;i<channels.size();i++) {
    sum += channels[i]->Alpha();
    if (sum>rn) {
      channels[i]->GeneratePoint(p);
      break;
    }
  }  
}

void Multi_Channel::GenerateWeight(double sprime,double y,int mode) {
  m_weight = 0.;
  for (short int i=0; i<channels.size(); ++i) {
    if (channels[i]->Alpha() > 0.) {
      channels[i]->GenerateWeight(sprime,y,mode);
      if (!(channels[i]->Weight()>0) && 
	  !(channels[i]->Weight()<0) && (channels[i]->Weight()!=0)) {
	AORGTOOLS::msg.Error()<<"Channel "<<i<<" produces a nan!"<<endl;
      }
      if (channels[i]->Weight()!=0.) 
	m_weight += channels[i]->Alpha()/channels[i]->Weight();
    }
  }
  if (!AMATOOLS::IsZero(m_weight)) m_weight = 1./m_weight;
}

void Multi_Channel::GenerateWeight(int n,double sprime,double y,int mode) {
  if (channels[n]->Alpha() > 0.) {
    channels[n]->GenerateWeight(sprime,y,mode);
    if (channels[n]->Weight()==0.) m_weight = 0.; 
                              else m_weight = 1./channels[n]->Weight();
  }
  else m_weight = 0.;
}

void Multi_Channel::GeneratePoint(double & sprime,double & y,int mode) {
  for(short int i=0;i<channels.size();i++) channels[i]->SetWeight(0.);
  double disc = ran.Get();
  double sum  = 0;
  for (short int n=0;n<channels.size();n++) {
    sum += channels[n]->Alpha();
    if (sum>disc) {
      for (int i=0;i<2;i++) rans[i] = ran.Get();
      channels[n]->GeneratePoint(sprime,y,mode,rans);
      return;
    }
  }  
}

void Multi_Channel::GeneratePoint(int n,double & sprime,double & y,int mode) {
  for (int i=0;i<2;i++) rans[i] = ran.Get();
  channels[n]->GeneratePoint(sprime,y,mode,rans);
}


int Multi_Channel::CountResonances(int i,Flavour*& fl_res) 
{
  return channels[i]->CountResonances(fl_res); 
}

void Multi_Channel::ISRInfo(int i,int & type,double & mass,double & width) 
{
  channels[i]->ISRInfo(type,mass,width);
  return;
}


void Multi_Channel::Print() {
  AORGTOOLS::msg.Out()<<"----------------------------------------------"<<endl
		      <<"Multi_Channel with "<<channels.size()<<" channels."<<endl;
  for (int i=0;i<channels.size();i++) 
    AORGTOOLS::msg.Out()<<"  "<<channels[i]->Name()<<" : "<<channels[i]->Alpha()<<endl;
  AORGTOOLS::msg.Out()<<"----------------------------------------------"<<endl;
}                 

void Multi_Channel::WriteOut(std::string pID) {
  ofstream ofile;
  ofile.open(pID.c_str());

  ofile<<channels.size()<<" "<<name<<endl;
  ofile.precision(12);
  for (int i=0;i<channels.size();i++) 
    ofile<<channels[i]->Name()<<" "<<channels[i]->Alpha()<<endl;
  ofile.close();
}

bool Multi_Channel::ReadIn(std::string pID) {
  ifstream ifile;
  ifile.open(pID.c_str());

  int         _size;
  std::string _name;
  double      _alpha;
  ifile>>_size>>_name;
  if (( _size != channels.size()) || ( _name != name) ) {
    msg.Error()<<"Error in Multi_Channel::ReadIn("<<pID<<")"<<endl 
	       <<"  Multi_Channel file did not coincide with actual Multi_Channel: "<<endl
	       <<"  "<<_size<<" vs. "<<channels.size()<<" and "
	       <<"  "<<_name<<" vs. "<<name<<endl;
    return 0;
  }

  double sum=0;
  for (int i=0;i<channels.size();i++) {
    ifile>>_name>>_alpha;
    sum+= _alpha;
    if (_name != channels[i]->Name()) {
      msg.Error()<<"Error in Multi_Channel::ReadIn("<<pID<<")"<<endl 
		 <<"  name of Single_Channel not consistent ("<<i<<")"<<endl
		 <<"  "<<_name<<" vs. "<<channels[i]->Name()<<endl;
      return 0;
    }
    channels[i]->SetAlpha(_alpha);
  }
  for (int i=0;i<channels.size();i++) {
    double a=channels[i]->Alpha();
    channels[i]->SetAlpha(a/sum);
  }

  ifile.close();
  return 1;
}
