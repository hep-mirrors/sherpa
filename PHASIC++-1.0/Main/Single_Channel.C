#include "Single_Channel.H"
#include "Random.H"
#include "MathTools.H"

using namespace PHASIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Single_Channel::Single_Channel(int _nin,int _nout,Flavour * fl) 
{ 
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout+1];
  for (short int i=0;i<nin+nout;i++) ms[i] = AMATOOLS::sqr(fl[i].Mass());

  if (nin == 1) rannum = 2 + 3*(nout-2);
  if (nin == 2) rannum = 1 + 2 + 3*(nout-2);
  rans  = new double[rannum];
}

Single_Channel::Single_Channel(Single_Channel * old)
{
  msg.Debugging()<<"New copy for Single_Channel : "<<old->name<<endl;
  nin    = old->nin;
  nout   = old->nout;
  rannum = old->rannum;
  rans   = new double[rannum];

  name   = old->name;

  ms     = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = old->ms[i];
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
  if (!AMATOOLS::IsZero(Value)) n_contrib++;
  n_points++;
  result  += Value;
  result2 += Value*Value;
};


void Single_Channel::GeneratePoint(Vec4D* p,Cut_Data * cuts)
{
  for (short int i=1;i<rannum;i++) rans[i] = ran.Get();
  //  cout<<"rannum="<<rannum<<endl;
  GeneratePoint(p,cuts,rans);
}


void Single_Channel::GeneratePoint(Vec4D * p)
{
  msg.Debugging()<<"In Single_Channel::GeneratePoint(Vec4d) "<<name<<endl;
  for (short int i=1;i<rannum;i++) rans[i] = ran.Get();
  GeneratePoint(p,rans);
}


