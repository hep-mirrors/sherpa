#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C4_5 : public Single_Channel {
    double m_amct,m_alpha,m_ctmax,m_ctmin;
    Info_Key m_kI_3_4,m_kTC_0__1_5_34_2,m_kTC_0__1__234_5,m_kZS_180;
    Vegas* p_vegas;
  public:
    C4_5(int,int,Flavour*,Integration_Info * const);
    ~C4_5();
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

extern "C" Single_Channel * Getter_C4_5(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C4_5(nin,nout,fl,info);
}

void C4_5::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p2345=p[0]+p[1];
  double s2345_max = p2345.Abs2();
  double s5 = ms[5];
  double s234_max = sqr(sqrt(s2345_max)-sqrt(ms[5]));
  double s2 = ms[2];
  double s25_min = cuts->Getscut(std::string("25"));
  double s34_max = sqr(sqrt(s2345_max)-sqrt(s25_min));
  double s3 = ms[3];
  double s4 = ms[4];
  double s34_min = cuts->Getscut(std::string("34"));
  Flavour fl34 = Flavour((kf_code)(24));
  Vec4D  p34;
  double s34 = CE.MassivePropMomenta(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,ran[0]);
  double s234_min = cuts->Getscut(std::string("234"));
  s234_min = Max(s234_min,sqr(sqrt(s34)+sqrt(s2)));
  Vec4D  p234;
  double s234 = CE.ThresholdMomenta(1.5,2.*sqrt(s234_min),s234_min,s234_max,ran[1]);
  m_ctmax = cuts->cosmax[1][5];
  m_ctmin = cuts->cosmin[1][5];
  CE.TChannelMomenta(p[0],p[1],p234,p[5],s234,s5,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[2],ran[3]);
  Vec4D  p1_5 = p[1]-p[5];
  CE.TChannelMomenta(p[0],p1_5,p34,p[2],s34,s2,0.,m_alpha,1.,-1.,m_amct,0,ran[4],ran[5]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[6],ran[7]);
}

void C4_5::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p2345=p[0]+p[1];
  double s2345_max = p2345.Abs2();
  double s5 = ms[5];
  double s234_max = sqr(sqrt(s2345_max)-sqrt(ms[5]));
  double s2 = ms[2];
  double s25_min = cuts->Getscut(std::string("25"));
  double s34_max = sqr(sqrt(s2345_max)-sqrt(s25_min));
  double s3 = ms[3];
  double s4 = ms[4];
  double s34_min = cuts->Getscut(std::string("34"));
  Flavour fl34 = Flavour((kf_code)(24));
  Vec4D  p34 = p[3]+p[4];
  double s34 = dabs(p34.Abs2());
  wt *= CE.MassivePropWeight(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,s34,rans[0]);
  double s234_min = cuts->Getscut(std::string("234"));
  s234_min = Max(s234_min,sqr(sqrt(s34)+sqrt(s2)));
  Vec4D  p234 = p[2]+p[3]+p[4];
  double s234 = dabs(p234.Abs2());
  wt *= CE.ThresholdWeight(1.5,2.*sqrt(s234_min),s234_min,s234_max,s234,rans[1]);
  m_ctmax = cuts->cosmax[1][5];
  m_ctmin = cuts->cosmin[1][5];
  if (m_kTC_0__1__234_5.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kTC_0__1__234_5<<CE.TChannelWeight(p[0],p[1],p234,p[5],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,m_kTC_0__1__234_5[0],m_kTC_0__1__234_5[1]);
  wt *= m_kTC_0__1__234_5.Weight();

  rans[2]= m_kTC_0__1__234_5[0];
  rans[3]= m_kTC_0__1__234_5[1];
  Vec4D  p1_5 = p[1]-p[5];
  if (m_kTC_0__1_5_34_2.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kTC_0__1_5_34_2<<CE.TChannelWeight(p[0],p1_5,p34,p[2],0.,m_alpha,1.,-1.,m_amct,0,m_kTC_0__1_5_34_2[0],m_kTC_0__1_5_34_2[1]);
  wt *= m_kTC_0__1_5_34_2.Weight();

  rans[4]= m_kTC_0__1_5_34_2[0];
  rans[5]= m_kTC_0__1_5_34_2[1];
  if (m_kI_3_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_3_4<<CE.Isotropic2Weight(p[3],p[4],m_kI_3_4[0],m_kI_3_4[1]);
  wt *= m_kI_3_4.Weight();

  rans[6]= m_kI_3_4[0];
  rans[7]= m_kI_3_4[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*3.-4.);

  weight = wt;
}

C4_5::C4_5(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C4_5");
  rannum = 8;
  rans  = new double[rannum];
  m_amct  = 1.;
  m_alpha = .9;
  m_ctmax = 1.;
  m_ctmin = -1.;
  m_kI_3_4.Assign(std::string("I_3_4"),2,0,info);
  m_kTC_0__1_5_34_2.Assign(std::string("TC_0__1_5_34_2"),2,0,info);
  m_kTC_0__1__234_5.Assign(std::string("TC_0__1__234_5"),2,0,info);
  m_kZS_180.Assign(std::string("ZS_180"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

C4_5::~C4_5()
{
  delete p_vegas;
}

void C4_5::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 180.94275;
  width = 0.;
}

void C4_5::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C4_5::ChID()
{
  return std::string("CG2$I_3_4$MP24_34$MTH_234$TC_0__1_5_34_2$TC_0__1__234_5$ZS_180$");
}
