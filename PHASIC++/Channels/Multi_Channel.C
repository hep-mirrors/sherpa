#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/My_File.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Multi_Channel::Multi_Channel(string _name) :
  name(StringReplace(_name, " ", "")),
  s1(NULL), m_readin(false), m_weight(1.0),
  n_points(0), n_contrib(0),
  mn_points(0), mn_contrib(0),
  m_lastdice(-1),
  m_otype(0)
{ }

Multi_Channel::~Multi_Channel()
{
  DropAllChannels();
  if (s1) { delete[] s1; s1 = NULL; }
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
    return;
  }
  if (s1!=NULL) delete[] s1;
  s1 =  new double[channels.size()];
  if (!m_readin) {
    s1xmin     = 1.e32;
    n_points   = 0;
    n_contrib  = 0;
    mn_points = mn_contrib = 0;
  }
  msg_Tracking()<<"Channels for "<<name<<" ("<<this<<")\n"
		<<"----------------- "<<n_points<<" --------------------"<<endl;
  for(size_t i=0;i<channels.size();i++) {
    if (!m_readin) channels[i]->Reset(1./channels.size());
    msg_Tracking()<<" "<<i<<" : "<<channels[i]->Name()<<"  : "<<channels[i]->Alpha()<<endl;
  }
  msg_Tracking()<<"----------------- "<<n_points<<" --------------------"<<endl;
  m_readin=false;
}

void Multi_Channel::MPISync()
{
#ifdef USING__MPI
  int size=mpi->Size();
  if (size>1) {
    int cn=2*channels.size()+2;
    double *values = new double[cn];
    for (size_t i=0;i<channels.size();++i) {
      values[i]=channels[i]->MRes1();
      values[channels.size()+i]=channels[i]->MRes2();
    }
    values[cn-2]=mn_points;
    values[cn-1]=mn_contrib;
    mpi->Allreduce(values,cn,MPI_DOUBLE,MPI_SUM);
    for (size_t i=0;i<channels.size();++i) {
      channels[i]->SetMPIVars(values[i],
			      values[channels.size()+i]);
    }
    mn_points=values[cn-2];
    mn_contrib=values[cn-1];
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
  msg_Tracking()<<"Optimize Multi_Channel : "<<name<<endl;

  size_t i;

  // calculate aptot = sum_i alpha_i
  double aptot = 0.;
  for (i=0;i<channels.size();i++) {
    s1[i]  = channels[i]->Res1()/n_points;
    aptot += channels[i]->Alpha()*sqrt(s1[i]);
  }

  // calculate s1x = max_i |aptot - sqrt(s1_i)|
  // update alpha_i -> alpha_i * sqrt(s1_i) / aptot
  //                     = alpha_i * sqrt(W_i(alpha_i))
  // where the last expression is given in the notation of hep-ph/9405257
  double s1x = 0.;
  for (i=0;i<channels.size();i++) {
    if (channels[i]->Alpha()>0.) {
      if (dabs(aptot-sqrt(s1[i]))>s1x) s1x = dabs(aptot-sqrt(s1[i]));
      if (channels.size()>1) {
        channels[i]->SetAlpha(channels[i]->Alpha() * sqrt(s1[i])/aptot);
      }
    }
  }

  // normalise alpha values to a partition of unity
  double norm = 0;
  for (i=0;i<channels.size();i++) norm += channels[i]->Alpha();
  for (i=0;i<channels.size();i++) {
    channels[i]->SetAlpha(channels[i]->Alpha() / norm);
  }

  // optimise individual channels ...
  for (i=0;i<channels.size();i++) channels[i]->Optimize();

  // save current alpha values if we have improved
  if (s1x<s1xmin) {
    s1xmin = s1x;
    for (i=0;i<channels.size();i++) {
      channels[i]->SetAlphaSave(channels[i]->Alpha());
    }
  }

  // reset channel weights
  for(i=0;i<channels.size();i++) channels[i]->ResetOpt();

  msg_Tracking()<<"New weights for : "<<name<<endl
		<<"----------------- "<<n_points<<" ----------------"<<endl;
  for (i=0;i<channels.size();i++) {
    if (channels[i]->Alpha() > 0) {
      msg_Tracking()<<i<<" channel "<<channels[i]->Name()<<" : "
		    <<channels[i]->Alpha()<<" -> "<<channels[i]->AlphaSave()<<endl;
    }
  }
  msg_Tracking()<<"S1X: "<<s1x<<" -> "<<s1xmin<<endl
 		<<"n,n_contrib : "<<n_points<<", "<<n_contrib<<endl
		<<"-----------------------------------------------"<<endl;
}

void Multi_Channel::EndOptimize(double error)
{
  size_t i;

  // use last best set of alpha values
  for (i=0;i<channels.size();i++) {
    channels[i]->SetAlpha(channels[i]->AlphaSave());
  }

  // normalise alpha values to a partition of unity
  double norm = 0;
  for (i=0;i<channels.size();i++) norm += channels[i]->Alpha();
  for (i=0;i<channels.size();i++) {
    channels[i]->SetAlpha(channels[i]->Alpha() / norm);
  }

  // tell channels to end optimising
  for (i=0;i<channels.size();i++) channels[i]->EndOptimize();

  msg_Tracking()<<"Best weights:-------------------------------"<<endl;
  for (i=0;i<channels.size();i++) {
    if (channels[i]->Alpha() > 0) {
      msg_Tracking()<<i<<" channel "<<channels[i]->Name()
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
  // update number of points
#ifdef USING__MPI
  if (value!=0.) mn_contrib++;
  mn_points++;
#else
  if (value!=0.) n_contrib++;
  n_points++;
#endif

  // update weights of all channels
  double var;
  for (size_t i=0;i<channels.size();i++) {
    if (value!=0.) {
      if (channels[i]->Weight()!=0) {
	var = sqr(value)*m_weight/channels[i]->Weight();
      } else {
        var = 0.;
      }
#ifdef USING__MPI
      channels[i]->AddMPIVars(var,sqr(var));
#else
      channels[i]->SetRes1(channels[i]->Res1() + var);
      channels[i]->SetRes2(channels[i]->Res2() + sqr(var));
#endif
    }
  }

  // add point to last selected channel
  if (m_lastdice>=0) Channel(m_lastdice)->AddPoint(value);
}

void Multi_Channel::GenerateWeight(Vec4D * p,Cut_Data * cuts)
{
  if (channels.empty()) return;
  Vec4D_Vector pp(p,&p[nin+nout]);
  m_weight = 0.;
  if (channels.size()==1) {
    channels[0]->GenerateWeight(&pp.front(),cuts);
    if (channels[0]->Weight()!=0) m_weight = channels[0]->Weight();
    return;
  }
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
  if (m_erans.size()) msg_Debugging()<<METHOD<<"(): Generating variables\n";
  for (std::map<std::string,double>::iterator
	 it(m_erans.begin());it!=m_erans.end();++it) {
    it->second=ran->Get();
    msg_Debugging()<<"  "<<it->first<<" -> "<<it->second<<"\n";
  }
  if (channels.empty()) {
    if (nin>1) p[2]=p[0]+p[1];
    else p[1]=p[0];
    return;
  }
  Poincare cms(p[0]+p[1]);
  for(size_t i=0;i<channels.size();i++) channels[i]->SetWeight(0.);
  if(channels.size()==1) {
    channels[0]->GeneratePoint(p,cuts);
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
      m_lastdice = i;
      break;
    }
  }
}

void Multi_Channel::GeneratePoint() {
  if (m_erans.size()) msg_Debugging()<<METHOD<<"(): Generating variables\n";
  for (std::map<std::string,double>::iterator
	 it(m_erans.begin());it!=m_erans.end();++it) {
    it->second=ran->Get();
    msg_Debugging()<<"  "<<it->first<<" -> "<<it->second<<"\n";
  }
  for(size_t i=0;i<channels.size();++i) channels[i]->SetWeight(0.);
  double disc=ran->Get();
  double sum=0.;
  for (size_t i=0;i<2;++i) rans[i]=ran->Get();
  for (size_t n=0;n<channels.size();++n) {
    sum += channels[n]->Alpha();
    if (sum>disc) {
      channels[n]->GeneratePoint(rans);
      m_lastdice = n;
      return;
    }
  }
  if (IsEqual(sum,disc)) {
    channels.back()->GeneratePoint(rans);
    m_lastdice = channels.size()-1;
    return;
  }
  msg_Error()<<"Multi_Channel::GeneratePoint("<<name<<"): IS case ("<<this
	     <<") No channel selected. \n"
	     <<"   disc = "<<disc<<", sum = "<<sum<<std::endl;
  Abort();
}

void Multi_Channel::GenerateWeight()
{
  if (channels.size()==1) {
    channels[0]->GenerateWeight();
    if (channels[0]->Weight()!=0) m_weight = channels[0]->Weight();
    return;
  }
  m_weight = 0.;
  for (size_t i=0;i<channels.size();++i) {
    if (channels[i]->Alpha()>0.) {
      channels[i]->GenerateWeight();
      if (!(channels[i]->Weight()>0)&&
	  !(channels[i]->Weight()<0)&&(channels[i]->Weight()!=0)) {
	msg_Error()<<"Multi_Channel::GenerateWeight(): ("<<this->name
		   <<"): Channel "<<i<<" ("<<channels[i]->Name()<<") produces a nan!"<<endl;
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
  My_Out_File ofile(pID);
  ofile.Open();
  ofile->precision(12);
  *ofile<<channels.size()<<" "<<name<<" "<<n_points<<" "<<n_contrib<<" "<<s1xmin<<endl;
//        <<m_result<<" "<<m_result2<<" "<<s1xmin<<" "
//        <<m_sresult<<" "<<m_sresult2<<" "<<m_ssigma2<<" "<<n_spoints<<" "<<m_optcnt<<endl;
  for (size_t i=0;i<channels.size();i++)
    *ofile<<channels[i]->Name()<<" "<<n_points<<" "
	 <<channels[i]->Alpha()<<" "<<channels[i]->AlphaSave()<<" "
	 <<0<<" "<<channels[i]->Res1()<<" "
	 <<channels[i]->Res2()<<std::endl;
  ofile.Close();
  for (size_t i=0;i<channels.size();i++) channels[i]->WriteOut(pID);
}

bool Multi_Channel::ReadIn(std::string pID) {
  My_In_File ifile(pID);
  if (!ifile.Open()) return false;
  size_t      size;
  std::string rname;
  long int    points;
  double      alpha, alphasave, weight, res1, res2;
  *ifile>>size>>rname;
  if (( size != channels.size()) || ( rname != name) ) {
    msg_Error()<<METHOD<<"(): Error reading in pID="<<pID<<endl
	       <<"  Multi_Channel file did not coincide with actual Multi_Channel: "<<endl
	       <<"  "<<size<<" vs. "<<channels.size()<<" and "
	       <<"  "<<rname<<" vs. "<<name<<endl;
    return 0;
  }
  m_readin=true;
  *ifile>>n_points>>n_contrib>>s1xmin;

  double sum=0;
  for (size_t i=0;i<channels.size();i++) {
    *ifile>>rname>>points>>alpha>>alphasave>>weight>>res1>>res2;
    sum+= alpha;
    if (rname != channels[i]->Name()) {
      msg_Error()<<METHOD<<"(): Error reading in pID="<<pID<<endl
		 <<"  name of Single_Channel not consistent ("<<i<<")"<<endl
		 <<"  "<<name<<" vs. "<<channels[i]->Name()<<endl;
      return 0;
      if (rname.substr(0,rname.length()-1)!=
          channels[i]->Name().substr(0,rname.length()-1)) {
	msg_Error()<<"   return 0."<<std::endl;
	return 0;
      }
    }
    channels[i]->SetAlpha(alpha);
    channels[i]->SetAlphaSave(alphasave);
    channels[i]->SetRes1(res1);
    channels[i]->SetRes2(res2);
  }
  ifile.Close();
  for (size_t i=0;i<channels.size();i++) channels[i]->ReadIn(pID);
  return 1;
}

std::string Multi_Channel::ChID(int n)
{
  return channels[n]->ChID();
}

bool Multi_Channel::Initialize()
{
  return true;
}
