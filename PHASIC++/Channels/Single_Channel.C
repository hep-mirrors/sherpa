#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/MathTools.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Single_Channel::Single_Channel() :
  name("no_name"),
  n_points(0),n_contrib(0),weight(0.),result(0.),result2(0.),
  res1(0.),res2(0.),res3(0.),alpha(0.),alpha_save(0.),
  nin(0),nout(0),ms(NULL),rannum(0),rans(NULL) 
{}

Single_Channel::Single_Channel(int _nin,int _nout,const Flavour * _fl) :
  name("no_name"),
  n_points(0),n_contrib(0),weight(0.),result(0.),result2(0.),
  res1(0.),res2(0.),res3(0.),alpha(0.),alpha_save(0.),
  nin(_nin),nout(_nout),ms(new double[nin+nout+1]),rannum(0),rans(NULL) 
{ 
  for (int i(0);i<nin+nout;i++) ms[i] = ATOOLS::sqr(_fl[i].Mass());
  rannum = 0;
  rans   = NULL;
  alpha=0.0;
  result=result2=0.0;
  //   if (nin == 1) rannum = 2 + 3*(nout-2);
  //   if (nin == 2) rannum = 1 + 2 + 3*(nout-2);
  //   rans  = new double[rannum];
}

Single_Channel::Single_Channel(Single_Channel * old) :
  name(old->name),
  n_points(0),n_contrib(0),weight(0.),result(0.),result2(0.),
  res1(0.),res2(0.),res3(0.),alpha(0.),alpha_save(0.),
  nin(old->nin),nout(old->nout),ms(new double[nin+nout]),
  rannum(old->rannum),rans(new double[rannum])
{
  for (int i=0;i<nin+nout;i++) ms[i] = old->ms[i];
  alpha=0.0;
  result=result2=0.0;
}

Single_Channel::~Single_Channel()
{
  if (ms) delete[] ms; 
  if (rans) delete[] rans; 
}

void Single_Channel::Reset(double value) {
  alpha    = alpha_save = value;
  weight   = 0.;
  res1     = res2       = res3 = 0.;
  result   = result2    = 0.0;
  n_points = n_contrib  = 0;
}

void Single_Channel::ResetOpt() {
  res1     = res2      = res3 = 0.;
  n_points = n_contrib = 0;
}

void Single_Channel::AddPoint(double Value) {
  //if (!ATOOLS::IsZero(Value)) n_contrib++;
  if (Value>0.) n_contrib++;
  n_points++;
  result  += Value;
  result2 += Value*Value;
}


void Single_Channel::GeneratePoint(Vec4D* p,Cut_Data * cuts)
{
  for (int i=0;i<rannum;i++) rans[i] = ran.Get();
  GeneratePoint(p,cuts,rans);
}


void Single_Channel::GeneratePoint(Vec4D * p)
{
  for (int i=0;i<rannum;i++) rans[i] = ran.Get();
  GeneratePoint(p,rans);
}

void Single_Channel::GeneratePoint(ATOOLS::Vec4D *p,Cut_Data *cuts,double *rans) 
{
  msg_Error()<<"Single_Channel::GeneratePoint(Vec4D *p,Cut_Data *cuts,double *rans): "
		     <<"Virtual Method called !"<<std::endl;
}

void Single_Channel::GenerateWeight(ATOOLS::Vec4D *p,Cut_Data *cuts) 
{
  msg_Error()<<"Single_Channel::GenerateWeight(Vec4D *p,Cut_Data *cuts): "
		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GeneratePoint(ATOOLS::Vec4D *p,double *rans) 
{ 
  msg_Error()<<"Single_Channel::GeneratePoint(Vec4D *p,double *rans).: "
		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GenerateWeight(ATOOLS::Vec4D *p) 
{
  msg_Error()<<"Single_Channel::GenerateWeight(Vec4D *p): "
		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GeneratePoint(double &sprime,double &y,int mode,double *rans) 
{
  msg_Error()<<"Single_Channel::GeneratePoint(double &sprime,double &y,int mode,double *rans): "
		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GeneratePoint(double &sprime,double &y,int mode) 
{
  msg_Error()<<"Single_Channel::GeneratePoint(double &sprime,double &y,int mode): "
		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GenerateWeight(double sprime,double y,int mode) 
{
  msg_Error()<<"Single_Channel::GenerateWeight(double sprime,double y,int mode): "
		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode) 
{
  msg_Error()<<"Single_Channel::GeneratePoint("<<mode<<"): "
		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const int mode) 
{
  msg_Error()<<"Single_Channel::GeneratePoint("<<mode<<"): "
 		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GenerateWeight(const int mode) 
{
  msg_Error()<<"Single_Channel::GenerateWeight("<<mode<<"): "
		     <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::SetRange(double *_sprimerange,double *_yrange) 
{
  for (int i=0;i<2;++i) {
    sprimerange[i]=_sprimerange[i];
    yrange[i]=_yrange[i];
  }
  sprimerange[2]=_sprimerange[2];
  m_Q=sqrt(sprimerange[2]);
}

void Single_Channel::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
{
  msg_Error()<<"Single_Channel::CalculateLimits(..): "
 		     <<"Virtual method called!"<<std::endl;
}

void Single_Channel::CalculateLimits() 
{
  msg_Error()<<"Single_Channel::CalculateLimits(): "
 		     <<"Virtual method called!"<<std::endl;
}

void Single_Channel::GetRange() 
{
  msg_Debugging()<<"  sprime : "<<sprimerange[0]<<" "<<sprimerange[1]<<" / "<<sprimerange[2]<<" / "
			 <<"  y : "<<yrange[0]<<" ... "<<yrange[1]<<std::endl;
}

void Single_Channel::ISRInfo(int &type,double &mass,double &width) 
{
  type=0;
  mass=width=0.0;
}


int Single_Channel::ChNumber() 
{
  msg_Error()<<"Method : Single_Channel::ChNumber()"<<std::endl;
  return 0;
}

void Single_Channel::SetChNumber(int) 
{
  msg_Error()<<"Method : Single_Channel::SetChNumber()"<<std::endl;
}

std::string Single_Channel::ChID() 
{ 
  msg_Error()<<"Virtual Method : Single_Channel::ChID()"<<std::endl;
  return std::string(""); 
}

size_t Single_Channel::Dimension() const 
{ 
  return rannum;
}

std::string      Single_Channel::Name()      { return name; }

int Single_Channel::N()         { return n_points; }
int Single_Channel::ValidN()    { return n_contrib; }
int Single_Channel::Nin()       { return nin; }
int Single_Channel::Nout()      { return nout; }

double Single_Channel::Res1()      { return res1; }
double Single_Channel::Res2()      { return res2; }
double Single_Channel::Res3()      { return res3; }
double Single_Channel::Weight()    { return weight; }
double Single_Channel::Alpha()     { return alpha; }
double Single_Channel::AlphaSave() { return alpha_save; }
double Single_Channel::Result()    { return result; }

void Single_Channel::IncrementN()                { n_points++; }
void Single_Channel::SetRes1(double _r)          { res1       = _r; }
void Single_Channel::SetRes2(double _r)          { res2       = _r; }
void Single_Channel::SetRes3(double _r)          { res3       = _r; }
void Single_Channel::SetN(long int _n)           { n_points   = _n; }
void Single_Channel::SetName(std::string _name)  { name       = _name; }
void Single_Channel::SetWeight(double _weight)   { weight     = _weight; }
void Single_Channel::SetAlpha(double _alpha)     { alpha      = _alpha; }
void Single_Channel::SetAlphaSave(double _alpha) { alpha_save = _alpha; }

void Single_Channel::Optimize() {}
void Single_Channel::EndOptimize() {}
void Single_Channel::WriteOut(std::string) {}
void Single_Channel::ReadIn(std::string) {}

int  Single_Channel::OType() { return 0; }

void Single_Channel::ISRInfo
(std::vector<int> &ts,std::vector<double> &ms,std::vector<double> &ws) const
{
}
