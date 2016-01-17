#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C3_5 : public Single_Channel {
    double m_amct,m_alpha,m_ctmax,m_ctmin;
    Info_Key m_kTC_0__1_3_4_2,m_kTC_0__1__24_3,m_kZS_0;
    Vegas* p_vegas;
  public:
    C3_5(int,int,Flavour*,Integration_Info * const);
    ~C3_5();
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

extern "C" Single_Channel * Getter_C3_5(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C3_5(nin,nout,fl,info);
}

void C3_5::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s3 = ms[3];
  double s24_max = sqr(sqrt(s234_max)-sqrt(ms[3]));
  double s2 = ms[2];
  double s4 = ms[4];
  double s24_min = cuts->Getscut(std::string("24"));
  Vec4D  p24;
  double s24 = CE.MasslessPropMomenta(.5,s24_min,s24_max,ran[0]);
  m_ctmax = cuts->cosmax[1][3];
  m_ctmin = cuts->cosmin[1][3];
  CE.TChannelMomenta(p[0],p[1],p24,p[3],s24,s3,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[1],ran[2]);
  Vec4D  p1_3 = p[1]-p[3];
  CE.TChannelMomenta(p[0],p1_3,p[4],p[2],s4,s2,0.,m_alpha,1.,-1.,m_amct,0,ran[3],ran[4]);
}

void C3_5::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s24_max = sqr(sqrt(s234_max)-sqrt(ms[3]));
  double s24_min = cuts->Getscut(std::string("24"));
  Vec4D  p24 = p[2]+p[4];
  double s24 = dabs(p24.Abs2());
  wt *= CE.MasslessPropWeight(.5,s24_min,s24_max,s24,rans[0]);
  m_ctmax = cuts->cosmax[1][3];
  m_ctmin = cuts->cosmin[1][3];
  if (m_kTC_0__1__24_3.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kTC_0__1__24_3<<CE.TChannelWeight(p[0],p[1],p24,p[3],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,m_kTC_0__1__24_3[0],m_kTC_0__1__24_3[1]);
  wt *= m_kTC_0__1__24_3.Weight();

  rans[1]= m_kTC_0__1__24_3[0];
  rans[2]= m_kTC_0__1__24_3[1];
  Vec4D  p1_3 = p[1]-p[3];
  if (m_kTC_0__1_3_4_2.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kTC_0__1_3_4_2<<CE.TChannelWeight(p[0],p1_3,p[4],p[2],0.,m_alpha,1.,-1.,m_amct,0,m_kTC_0__1_3_4_2[0],m_kTC_0__1_3_4_2[1]);
  wt *= m_kTC_0__1_3_4_2.Weight();

  rans[3]= m_kTC_0__1_3_4_2[0];
  rans[4]= m_kTC_0__1_3_4_2[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,3*3.-4.);

  weight = wt;
}

C3_5::C3_5(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C3_5");
  rannum = 5;
  rans  = new double[rannum];
  m_amct  = 1.0+ToType<double>(rpa->gen.Variable("AMEGIC_CHANNEL_EPSILON"));
  m_alpha = ToType<double>(rpa->gen.Variable("AMEGIC_CHANNEL_ALPHA"));
  m_ctmax = 1.;
  m_ctmin = -1.;
  m_kTC_0__1_3_4_2.Assign(std::string("TC_0__1_3_4_2"),2,0,info);
  m_kTC_0__1__24_3.Assign(std::string("TC_0__1__24_3"),2,0,info);
  m_kZS_0.Assign(std::string("ZS_0"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

C3_5::~C3_5()
{
  delete p_vegas;
}

void C3_5::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 0;
  width = 0.;
}

void C3_5::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C3_5::ChID()
{
  return std::string("CGND$MTH_24$TC_0__1_3_4_2$TC_0__1__24_3$ZS_0$");
}
