#include "Single_Channel.H"
#include "Random.H"
#include "MathTools.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Single_Channel::Single_Channel(int _nin,int _nout,Flavour * fl) 
{ 
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout+1];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());

  if (nin == 1) rannum = 2 + 3*(nout-2);
  if (nin == 2) rannum = 1 + 2 + 3*(nout-2);
  rans  = new double[rannum];
}

Single_Channel::Single_Channel(Single_Channel * old)
{
  nin    = old->nin;
  nout   = old->nout;
  rannum = old->rannum;
  rans   = new double[rannum];

  name   = old->name;

  ms     = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = old->ms[i];
}

Single_Channel::~Single_Channel()
{
  if (ms)   { delete[] ms; }
  if (rans) { delete[] rans; }
}

void Single_Channel::Reset(double value) {
  alpha    = alpha_save = value;
  weight   = 0.;
  res1     = res2       = res3 = 0.;
  n_points = n_contrib  = 0;
}

void Single_Channel::ResetOpt() {
  res1     = res2      = res3 = 0.;
  n_points = n_contrib = 0;
};

void Single_Channel::AddPoint(double Value) {
  if (!ATOOLS::IsZero(Value)) n_contrib++;
  n_points++;
  result  += Value;
  result2 += Value*Value;
};


void Single_Channel::GeneratePoint(Vec4D* p,Cut_Data * cuts)
{
  for (short int i=1;i<rannum;i++) rans[i] = ran.Get();
  GeneratePoint(p,cuts,rans);
}


void Single_Channel::GeneratePoint(Vec4D * p)
{
  for (short int i=1;i<rannum;i++) rans[i] = ran.Get();
  GeneratePoint(p,rans);
}

void Single_Channel::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode) 
{
  ATOOLS::msg.Error()<<"Single_Channel::GeneratePoint("<<mode<<"): "
                   <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const int mode) 
{
  ATOOLS::msg.Error()<<"Single_Channel::GeneratePoint("<<mode<<"): "
                   <<"Virtual Method called !"<<std::endl; 
}

void Single_Channel::GenerateWeight(const int mode) 
{
  ATOOLS::msg.Error()<<"Single_Channel::GenerateWeight("<<mode<<"): "
                   <<"Virtual Method called !"<<std::endl; 
}

bool Single_Channel::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
{
  ATOOLS::msg.Error()<<"Single_Channel::CalculateLimits(): "
                   <<"Virtual method called!"<<std::endl;
  return false;
}


void Single_Channel::GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *,double *) 
{
  ATOOLS::msg.Error()<<"Virtual Method : "
		     <<"Single_Channel::GeneratePoint(Vec4D,Cut_Data,double)."<<std::endl; 
}

void Single_Channel::GenerateWeight(ATOOLS::Vec4D *,ATOOLS::Cut_Data *) 
{
  ATOOLS::msg.Error()<<"Virtual Method : "
		     <<"Single_Channel::GenerateWeight(Vec4D,Cut_Data)"<<std::endl; 
}

void Single_Channel::GeneratePoint(ATOOLS::Vec4D *,double *) 
{ 
  ATOOLS::msg.Error()<<"Virtual Method : "
		     <<"Single_Channel::GeneratePoint(Vec4D,double)."<<std::endl; 
}

void Single_Channel::GenerateWeight(ATOOLS::Vec4D *) 
{
  ATOOLS::msg.Error()<<"Virtual Method : "
		     <<"Single_Channel::GenerateWeight(Vec4D)"<<std::endl; 
}

void Single_Channel::GeneratePoint(double &,double &,int,double *) 
{
  ATOOLS::msg.Error()<<"Virtual Method : "
		     <<"Single_Channel::GeneratePoint(double,double,double)"<<std::endl; 
}

void Single_Channel::GeneratePoint(double &,double &,int) 
{
  ATOOLS::msg.Error()<<"Virtual Method : "
		     <<"Single_Channel::GeneratePoint(double,double)"<<std::endl; 
}

void Single_Channel::GenerateWeight(double,double,int) 
{
  ATOOLS::msg.Error()<<"Virtual Method : "
		     <<"Single_Channel::GenerateWeight(double,double)"<<std::endl; 
}

void Single_Channel::SetRange(double * _sprimerange,double * _yrange) 
{
  for (int i=0;i<2;i++) {
    sprimerange[i] = _sprimerange[i];
    yrange[i]      = _yrange[i];
  }
  sprimerange[2] =  _sprimerange[2];
}

void Single_Channel::GetRange() 
{
  ATOOLS::msg.Debugging()<<"  sprime : "<<sprimerange[0]<<" "<<sprimerange[1]<<" / "<<sprimerange[2]<<" / "
			 <<"  y : "<<yrange[0]<<" ... "<<yrange[1]<<std::endl;
}

void Single_Channel::ISRInfo(int &,double &,double &) 
{
  ATOOLS::msg.Error()<<"Virtual Method : Single_Channel::ISRInfo()"<<std::endl;
}

int Single_Channel::CountResonances(ATOOLS::Flavour*&) 
{ 
  ATOOLS::msg.Error()<<"Virtual Method : Single_Channel::CountResonances()"<<std::endl; 
  return 0;
}

int Single_Channel::ChNumber() 
{
  ATOOLS::msg.Error()<<"Virtual Method : Single_Channel::ChNumber()"<<std::endl;
  return 0;
}

void Single_Channel::SetChNumber(int) 
{
  ATOOLS::msg.Error()<<"Virtual Method : Single_Channel::SetChNumber()"<<std::endl;
}

std::string Single_Channel::ChID() 
{ 
  ATOOLS::msg.Error()<<"Virtual Method : Single_Channel::ChID()"<<std::endl;
  return std::string(""); 
}
