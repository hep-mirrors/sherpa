#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C4_18 : public Single_Channel {
    Info_Key m_kI_2_345,m_kI_3_4,m_kI_5_34,m_kZS_180;
    Vegas* p_vegas;
  public:
    C4_18(int,int,Flavour*,Integration_Info * const);
    ~C4_18();
    void   GenerateWeight(Vec4D *,Cut_Data *);
    void   GeneratePoint(Vec4D *,Cut_Data *,double *);
    void   AddPoint(double);
    void   MPISync()                 { p_vegas->MPISync(); }
    void   Optimize()                { p_vegas->Optimize(); } 
    void   EndOptimize()             { p_vegas->EndOptimize(); } 
    void   WriteOut(std::string pId) { p_vegas->WriteOut(pId); } 
    void   ReadIn(std::string pId)   { p_vegas->ReadIn(pId); } 
    void   ISRInfo(int &,double &,double &);
    std::string ChID();
  };
}

extern "C" Single_Channel * Getter_C4_18(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C4_18(nin,nout,fl,info);
}

void C4_18::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p2345=p[0]+p[1];
  double s2345_max = p2345.Abs2();
  double s2 = ms[2];
  double s345_max = sqr(sqrt(s2345_max)-sqrt(ms[2]));
  double s25_min = cuts->Getscut(std::string("25"));
  double s34_max = sqr(sqrt(s2345_max)-sqrt(s25_min));
  double s3 = ms[3];
  double s4 = ms[4];
  double s34_min = cuts->Getscut(std::string("34"));
  Flavour fl34 = Flavour((kf_code)(24));
  Vec4D  p34;
  double s34 = CE.MassivePropMomenta(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,ran[0]);
  double s5 = ms[5];
  double s345_min = cuts->Getscut(std::string("345"));
  s345_min = Max(s345_min,sqr(sqrt(s5)+sqrt(s34)));
  Vec4D  p345;
  double s345 = CE.ThresholdMomenta(1.5,2.*sqrt(s345_min),s345_min,s345_max,ran[1]);
  CE.Isotropic2Momenta(p2345,s2,s345,p[2],p345,ran[2],ran[3]);
  CE.Isotropic2Momenta(p345,s5,s34,p[5],p34,ran[4],ran[5]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[6],ran[7]);
}

void C4_18::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p2345=p[0]+p[1];
  double s2345_max = p2345.Abs2();
  double s2 = ms[2];
  double s345_max = sqr(sqrt(s2345_max)-sqrt(ms[2]));
  double s25_min = cuts->Getscut(std::string("25"));
  double s34_max = sqr(sqrt(s2345_max)-sqrt(s25_min));
  double s3 = ms[3];
  double s4 = ms[4];
  double s34_min = cuts->Getscut(std::string("34"));
  Flavour fl34 = Flavour((kf_code)(24));
  Vec4D  p34 = p[3]+p[4];
  double s34 = dabs(p34.Abs2());
  wt *= CE.MassivePropWeight(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,s34,rans[0]);
  double s5 = ms[5];
  double s345_min = cuts->Getscut(std::string("345"));
  s345_min = Max(s345_min,sqr(sqrt(s5)+sqrt(s34)));
  Vec4D  p345 = p[3]+p[4]+p[5];
  double s345 = dabs(p345.Abs2());
  wt *= CE.ThresholdWeight(1.5,2.*sqrt(s345_min),s345_min,s345_max,s345,rans[1]);
  if (m_kI_2_345.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_2_345<<CE.Isotropic2Weight(p[2],p345,m_kI_2_345[0],m_kI_2_345[1]);
  wt *= m_kI_2_345.Weight();

  rans[2]= m_kI_2_345[0];
  rans[3]= m_kI_2_345[1];
  if (m_kI_5_34.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_5_34<<CE.Isotropic2Weight(p[5],p34,m_kI_5_34[0],m_kI_5_34[1]);
  wt *= m_kI_5_34.Weight();

  rans[4]= m_kI_5_34[0];
  rans[5]= m_kI_5_34[1];
  if (m_kI_3_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_3_4<<CE.Isotropic2Weight(p[3],p[4],m_kI_3_4[0],m_kI_3_4[1]);
  wt *= m_kI_3_4.Weight();

  rans[6]= m_kI_3_4[0];
  rans[7]= m_kI_3_4[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*3.-4.);

  weight = wt;
}

C4_18::C4_18(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C4_18");
  rannum = 8;
  rans  = new double[rannum];
  m_kI_2_345.Assign(std::string("I_2_345"),2,0,info);
  m_kI_3_4.Assign(std::string("I_3_4"),2,0,info);
  m_kI_5_34.Assign(std::string("I_5_34"),2,0,info);
  m_kZS_180.Assign(std::string("ZS_180"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

C4_18::~C4_18()
{
  delete p_vegas;
}

void C4_18::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 180.94275;
  width = 0.;
}

void C4_18::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C4_18::ChID()
{
  return std::string("CG2$I_2_345$I_3_4$I_5_34$MP24_34$MTH_345$ZS_180$");
}
