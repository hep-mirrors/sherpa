#include "Multi_Channel.H"
#include "Random.H"

using namespace PHASIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Multi_Channel::Multi_Channel(string _name) 
{
  string help;
  int    pos;
  for (;;) {
    pos  = _name.find(" ");
    if (pos==-1) break;
    help  = _name;
    _name = help.substr(0,pos) + help.substr(pos+1); 
  }
  name     = _name;
  s1       = s2        = 0;
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
  result     = 0.;
  result2    = 0.;

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
};        

/*
void Multi_Channel::MPIOptimize()
{
#ifdef _USE_MPI_
  int rank = MPI::COMM_WORLD.Get_Rank();
  int size = MPI::COMM_WORLD.Get_size();
  
  cout<<"Process "<<rank<<" in MPIOptimize()."<<endl;

  double* messageblock = new double[1+3*Chp.size()];
  double* alp = new double[Chp.size()];
  //tag 9 for communication
  if (rank==0) {
    int count = 0;
    while (count<size-1) {
      MPI::COMM_WORLD.Recv(messageblock, 1+3*Chp.size(), MPI::DOUBLE, MPI::ANY_SOURCE, 9);
      count++;
      cout<<" received Channel-Data from Knot "<<messageblock[0]<<endl;
      for (short int i=0;i<Chp.size();i++) {
	Chp[i]->Res1 += messageblock[1+Chp.size()+i];
	Chp[i]->Res2 += messageblock[1+2*Chp.size()+i];
	Chp[i]->Res3  = sqr(Chp[i]->Res1)-Chp[i]->Res2;
	Chp[i]->N    += int(messageblock[1+i]);
      }
    }
    cout<<"Master: "<<Chp[0]->N<<" data accumulated."<<endl; 
    Optimize();
    //broadcast alphai
    for (short int i=0;i<Chp.size();i++) alp[i] = Chp[i]->Alpha;
    
    for (int i=1;i<size;i++) MPI::COMM_WORLD.Isend(alp, Chp.size(), MPI::DOUBLE, i, 9);
  }
  else {
    messageblock[0] = rank;
    for (int i=0;i<Chp.size();i++) messageblock[1+i]       = Chp[i]->N;
    for (int i=0;i<Chp.size();i++) messageblock[1+Chp.size()+i]   = Chp[i]->Res1;
    for (int i=0;i<Chp.size();i++) messageblock[1+2*Chp.size()+i] = Chp[i]->Res2;
    MPI::COMM_WORLD.Send(messageblock, 1+3*Chp.size(), MPI::DOUBLE, 0, 9);   
    //Waiting for new alpha's
    MPI::COMM_WORLD.Recv(alp, Chp.size(), MPI::DOUBLE, 0, 9);
    cout<<"Slave "<<rank<<" received new alpha."<<endl;

    for (short int i=0;i<Chp.size();i++) {
      Chp[i]->Alpha = alp[i];
      Chp[i]->Reset_Opt();
      Chp[i]->Weight = 0.;
    }
  }
  delete[] alp;
  delete[] messageblock;
#else
  cout<<"MPIOptimize called in non-MPI session!!!"<<endl;
#endif

}
*/

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
  msg.Tracking()<<"result,result2,n,n_contrib : "<<result<<", ";
  msg.Tracking()<<result2<<", "<<n_points<<", "<<n_contrib<<endl;
  msg.Tracking()<<"-----------------------------------------------"<<endl;
}

void Multi_Channel::EndOptimize(double error)
{
  short int i;
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
  msg.Tracking()<<"result,result2,n,n_contrib : "<<result<<", ";
  msg.Tracking()<<result2<<", "<<n_points<<", "<<n_contrib<<endl;
  msg.Tracking()<<"-------------------------------------------"<<endl;
}

void Multi_Channel::AddPoint(double value)
{
  // msg.Debugging()<<"In Multi_Channel::AddPoint("<<value<<")"<<endl;
  if (!AMATOOLS::IsZero(value)) n_contrib++;

  n_points++;
  result  += value;
  result2 += value*value;
  
  double var;
  for (short int i=0;i<channels.size();i++) {
    if (value!=0.) {
      if (channels[i]->Weight()!=0) 
	var = sqr(value)*weight/channels[i]->Weight();
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
  double disc = (n_points*result2)/((n_points-1)*AMATOOLS::sqr(result)) - 1./(n_points-1);
  if (disc>0.) return result/n_points * sqrt(disc);
  disc = result2/(n_points*(n_points-1)) - AMATOOLS::sqr(result/n_points)/(n_points-1);
  if (disc>0.) return sqrt(disc);
  
  AORGTOOLS::msg.Error()<<"Variance yielded a NaN !"<<endl;
  AORGTOOLS::msg.Error()<<"   res,res2 = "<<result<<", "<<result2;
  AORGTOOLS::msg.Error()<<" after "<<n_points<<" points."<<endl; 
  
  return sqrt(-disc);
};


void Multi_Channel::GenerateWeight(int n,Vec4D* p,Cut_Data * cuts) 
{
  if (channels[n]->Alpha() > 0.) {
    channels[n]->GenerateWeight(p,cuts);
    if (channels[n]->Weight()==0.) weight = 0.; 
                              else weight = 1./channels[n]->Weight();
  }
  else weight = 0.;
}

void Multi_Channel::GenerateWeight(Vec4D * p,Cut_Data * cuts)
{
  weight = 0.;
  for (short int i=0; i<channels.size(); ++i) {
    if (channels[i]->Alpha() > 0.) {
      channels[i]->GenerateWeight(p,cuts);
      if (!(channels[i]->Weight()>0) && 
	  !(channels[i]->Weight()<0) && (channels[i]->Weight()!=0)) {
	msg.Error()<<"Channel "<<i<<" produces a nan!"<<endl;
      }
      if (channels[i]->Weight()!=0) 
	weight += channels[i]->Alpha()/channels[i]->Weight();
    }
  }
  if (weight!=0) weight = 1./weight;
}


void Multi_Channel::GeneratePoint(int n,Vec4D * p,Cut_Data * cuts,double * ran)
{
  channels[n]->GeneratePoint(p,cuts,ran);
}

void Multi_Channel::GeneratePoint(Vec4D* p,Cut_Data * cuts)
{
  for(short int i=0;i<channels.size();i++) channels[i]->SetWeight(0.);
  double rn  = ran.Get();
  double sum = 0;
  for (short int i=0;i<channels.size();i++) {
    sum += channels[i]->Alpha();
    if (sum>rn) {
      //      cout<<"Channel number "<<i<<endl;
      channels[i]->GeneratePoint(p,cuts);
      break;
    }
  }  
}

void Multi_Channel::GenerateWeight(int n,Vec4D* p) 
{
  if (channels[n]->Alpha() > 0.) {
    channels[n]->GenerateWeight(p);
    if (channels[n]->Weight()==0.) weight = 0.; 
                              else weight = 1./channels[n]->Weight();
  }
  else weight = 0.;
}


void Multi_Channel::GenerateWeight(Vec4D * p)
{
  weight = 0.;
  for (short int i=0; i<channels.size(); ++i) {
    if (channels[i]->Alpha() > 0.) {
      channels[i]->GenerateWeight(p);
      if (!(channels[i]->Weight()>0) && 
	  !(channels[i]->Weight()<0) && (channels[i]->Weight()!=0)) {
	msg.Error()<<"Channel "<<i<<" produces a nan!"<<endl;
      }
      if (channels[i]->Weight()!=0) 
	weight += channels[i]->Alpha()/channels[i]->Weight();
    }
  }
  if (!AMATOOLS::IsZero(weight)) weight = 1./weight;
}


void Multi_Channel::GeneratePoint(int n,Vec4D * p,double * rn)
{
  channels[n]->GeneratePoint(p,rn);
}

void Multi_Channel::GeneratePoint(Vec4D* p)
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
  weight = 0.;
  for (short int i=0; i<channels.size(); ++i) {
    if (channels[i]->Alpha() > 0.) {
      channels[i]->GenerateWeight(sprime,y,mode);
      if (!(channels[i]->Weight()>0) && 
	  !(channels[i]->Weight()<0) && (channels[i]->Weight()!=0)) {
	AORGTOOLS::msg.Error()<<"Channel "<<i<<" produces a nan!"<<endl;
      }
      if (channels[i]->Weight()!=0.) 
	weight += channels[i]->Alpha()/channels[i]->Weight();
    }
  }
  if (!AMATOOLS::IsZero(weight)) weight = 1./weight;
}

void Multi_Channel::GenerateWeight(int n,double sprime,double y,int mode) {
  if (channels[n]->Alpha() > 0.) {
    channels[n]->GenerateWeight(sprime,y,mode);
    if (channels[n]->Weight()==0.) weight = 0.; 
                              else weight = 1./channels[n]->Weight();
  }
  else weight = 0.;
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
  for (int i=0;i<channels.size();i++) {
    ifile>>_name>>_alpha;
    if (_name != channels[i]->Name()) {
      msg.Error()<<"Error in Multi_Channel::ReadIn("<<pID<<")"<<endl 
		 <<"  name of Single_Channel not consistent ("<<i<<")"<<endl
		 <<"  "<<_name<<" vs. "<<channels[i]->Name()<<endl;
      return 0;
    }
    channels[i]->SetAlpha(_alpha);
  }
  ifile.close();
  return 1;
}
