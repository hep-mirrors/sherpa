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
