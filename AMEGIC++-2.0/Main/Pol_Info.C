#include "Pol_Info.H"
#include "MathTools.H"
#include "Message.H"

using namespace AMEGIC;

Pol_Info::Pol_Info(const Pol_Info & p) 
{
  num=p.num;
  p_type=p.p_type;
  angle=p.angle;
  type=new int[num];
  factor=new double[num];
  for(int i=0;i<num;i++){
    type[i]=p.type[i];
    factor[i]=p.factor[i];
  }
}

Pol_Info& Pol_Info::operator=(const Pol_Info& p)
{
  if (this!=&p) {
    num    = p.num;
    p_type = p.p_type;
    angle  = p.angle;
    if(type)   delete[] type;
    if(factor) delete[] factor;
    type   = new int[num];
    factor = new double[num];
    for(int i=0;i<num;i++){
      type[i]   = p.type[i];
      factor[i] = p.factor[i];
    }
  }
  return *this;
}
Pol_Info::Pol_Info() { num=0; type=0; factor=0; p_type=' '; angle=0.;}
Pol_Info::Pol_Info(const APHYTOOLS::Flavour& fl)
{
  int dof = 1;
  if(fl.IsFermion())                                 { dof = 2;p_type='h';};
  if(fl.IsVector() &&  AMATOOLS::IsZero(fl.Mass()))  { dof = 2;p_type='c';}
  if(fl.IsVector() && !AMATOOLS::IsZero(fl.Mass()))  {

#ifdef Explicit_Pols
    dof=3;
#else
    dof=1;
#endif
    p_type='c';
  }
  Init((int)dof);
  int tf[3]  = {mt::p_m, mt::p_p, mt::p_l };
  for(int j=0;j<dof;j++){
    type[j]   = tf[j];
    factor[j] = 1.;
  }
}

Pol_Info::~Pol_Info(){if(type) delete[] type;if(factor)delete[] factor;}

void Pol_Info::Init(int i){num=i;type=new int[num];factor=new double[num];}



void Tensor_Struc::GetPolCombos(int num, std::vector<std::vector<int> >* pol, std::vector<int>* sign)
{
  pol->clear();
  sign->clear();
  std::vector<int> cc;cc.push_back(8);cc.push_back(8); 
  sign->push_back(1);
  switch(num){
  case mt::p_t1:
    cc[0]=mt::p_p;cc[1]=mt::p_p;
    pol->push_back(cc);
    break;
  case mt::p_t2:
    cc[0]=mt::p_p;cc[1]=mt::p_l;
    pol->push_back(cc);
    break;
  case mt::p_t3:
    cc[0]=mt::p_p;cc[1]=mt::p_m;
    pol->push_back(cc);
    cc[0]=mt::p_l;cc[1]=mt::p_l;
    pol->push_back(cc);
    sign->push_back(-1);
    break;
  case mt::p_t4:
    cc[0]=mt::p_m;cc[1]=mt::p_l;
    pol->push_back(cc);
    break;
  case mt::p_t5:
    cc[0]=mt::p_m;cc[1]=mt::p_m;
    pol->push_back(cc);
    break;
  default: AORGTOOLS::msg.Error()<<"Invalid tensor type: "<<num<<std::endl;
    abort();
  }
}

double Tensor_Struc::GetTfactor(int num)
{
  switch(num){
  case mt::p_t1:return 1.;
  case mt::p_t2:return 2.;
  case mt::p_t3:return 2./3.;
  case mt::p_t4:return 2.;
  case mt::p_t5:return 1.;
  }
  return 1.;
}
