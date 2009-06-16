#include "AMEGIC++/Main/Pol_Info.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;

std::ostream & AMEGIC::operator<<(std::ostream & s, Pol_Info & pi) 
{
  s<<" Pol_Info : "<<pi.pol_type<<std::endl;
  for (int i=0;i<pi.num;++i) {
    s<<pi.type[i]<<":"<<pi.factor[i]<<std::endl;
  }
  return s;
}

Pol_Info::Pol_Info(const Pol_Info & p) 
{
  num=p.num;
  pol_type=p.pol_type;
  angle=p.angle;
  if (num>0) {
    type=new int[num];
    factor=new double[num];
    for(int i=0;i<num;i++){
      type[i]=p.type[i];
      factor[i]=p.factor[i];
    }
  }
  else {
    type = 0;
    factor= 0;
  }
}

Pol_Info& Pol_Info::operator=(const Pol_Info& p)
{
  if (this!=&p) {
    num    = p.num;
    pol_type = p.pol_type;
    angle  = p.angle;
    if(type)   delete[] type;
    if(factor) delete[] factor;
    if (num>0) {
      type   = new int[num];
      factor = new double[num];
      for(int i=0;i<num;i++){
	type[i]   = p.type[i];
	factor[i] = p.factor[i];
      }
    }
    else {
      type = 0;
      factor= 0;
    }
  }
  return *this;
}
Pol_Info::Pol_Info() { num=0; type=0; factor=0; pol_type=' '; angle=0.;}
Pol_Info::Pol_Info(const ATOOLS::Flavour& fl)
{
  int dof = 1;
  pol_type='s';
  if(fl.IsFermion())                                 { dof = 2;pol_type='h';};
  if(fl.IsVector() &&  ATOOLS::IsZero(fl.Mass()))  { dof = 2;pol_type='c';}
  if(fl.IsVector() && !ATOOLS::IsZero(fl.Mass()))  {

#ifdef Explicit_Pols
    dof=3;
#else
    dof=1;
#endif
    pol_type='c';
  }
  if(fl.IsTensor()) dof=5;
  Init(dof);

  if(!fl.IsTensor()) {
    int tf[3]  = {mt::p_m, mt::p_p, mt::p_l };
    for(int j=0;j<dof;j++){
      type[j]   = tf[j];
      factor[j] = 1.;
    }
  } 
  else {
    type[0]=mt::p_t1;factor[0]=1.;
    type[1]=mt::p_t2;factor[1]=1.;
    type[2]=mt::p_t3;factor[2]=1.;
    type[3]=mt::p_t4;factor[3]=1.;
    type[4]=mt::p_t5;factor[4]=1.;
  }
}

Pol_Info::~Pol_Info(){if(type) delete[] type;if(factor)delete[] factor;}

void Pol_Info::Init(int i){num=i;type=new int[num];factor=new double[num];}

void Pol_Info::SetPol(char pol) 
{
  // fixes polarisation to + - or 0 of an unpolarised state:
  mt::momtype t1=mt::p_m, t2=mt::p_p;
  if (pol=='l') { 
    t1=mt::p_l0; 
    t2=mt::p_l1; 
  }
  mt::momtype t;
  switch (pol) {
    case '-': t = t1;     break;
    case '+': t = t2;     break;
    case '0': t = mt::p_l;break;
    default : t = t1;
  }
  int dof=num;
  if (!type) Init(1);
  type[0]=t;
  factor[0]=dof;
  num=1;
}

char Pol_Info::GetPol() {
  if (num==1) {
    if (pol_type=='s' || pol_type==' ') return 's';
    switch (type[0]) {
    case mt::p_m: return '-';
    case mt::p_p: return '+';
    case mt::p_l: return '0';
    case mt::p_l0: return 'x';
    case mt::p_l1: return 'y';

    case mt::p_t1: return 'a';
    case mt::p_t2: return 'b';
    case mt::p_t3: return 'c';
    case mt::p_t4: return 'd';
    case mt::p_t5: return 'e';
    default:
      return ' ';
    }
  }
  return ' ';
}

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
  default: 
    msg_Error()<<"ERROR in Tensor_Struc::GetPolCombos : "<<std::endl
		       <<"   Invalid tensor type: "<<num<<", abort the run."<<std::endl;
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
