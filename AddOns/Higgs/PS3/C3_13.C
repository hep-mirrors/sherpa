#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C3_13 : public Single_Channel {
    double m_alpha,m_ctmax,m_ctmin;
    Info_Key m_kI_3_4,m_kTC_0__1__34_2,m_kZS_0;
    Vegas* p_vegas;
  public:
    C3_13(int,int,Flavour*,Integration_Info * const);
    ~C3_13();
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

extern "C" Single_Channel * Getter_C3_13(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C3_13(nin,nout,fl,info);
}

void C3_13::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<m_rannum;i++) p_rans[i]=ran[i];
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s2 = p_ms[2];
  double s34_max = sqr(sqrt(s234_max)-sqrt(p_ms[2]));
  double s4 = p_ms[4];
  double s3 = p_ms[3];
  double s34_min = cuts->Getscut((1<<3)|(1<<4));
  Vec4D  p34;
  double s34 = CE.MasslessPropMomenta(.5,s34_min,s34_max,ran[0]);
  CE.TChannelMomenta(p[0],p[1],p34,p[2],s34,s2,0.,m_alpha,m_ctmax,m_ctmin,ran[1],ran[2]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[3],ran[4]);
}

void C3_13::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s34_max = sqr(sqrt(s234_max)-sqrt(p_ms[2]));
  double s34_min = cuts->Getscut((1<<3)|(1<<4));
  Vec4D  p34 = p[3]+p[4];
  double s34 = dabs(p34.Abs2());
  wt *= CE.MasslessPropWeight(.5,s34_min,s34_max,s34,p_rans[0]);
  if (m_kTC_0__1__34_2.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kTC_0__1__34_2<<CE.TChannelWeight(p[0],p[1],p34,p[2],0.,m_alpha,m_ctmax,m_ctmin,m_kTC_0__1__34_2[0],m_kTC_0__1__34_2[1]);
  wt *= m_kTC_0__1__34_2.Weight();

  p_rans[1]= m_kTC_0__1__34_2[0];
  p_rans[2]= m_kTC_0__1__34_2[1];
  if (m_kI_3_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_3_4<<CE.Isotropic2Weight(p[3],p[4],m_kI_3_4[0],m_kI_3_4[1]);
  wt *= m_kI_3_4.Weight();

  p_rans[3]= m_kI_3_4[0];
  p_rans[4]= m_kI_3_4[1];
  double vw = p_vegas->GenerateWeight(p_rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,3*3.-4.);

  m_weight = wt;
}

C3_13::C3_13(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  m_name = std::string("C3_13");
  m_rannum = 5;
  p_rans  = new double[m_rannum];
  Settings& s = Settings::GetMainSettings();
  m_alpha = s["SCHANNEL_ALPHA"].Get<double>();
  m_ctmax = 1.;
  m_ctmin = -1.;
  m_kI_3_4.Assign(std::string("I_3_4"),2,0,info);
  m_kTC_0__1__34_2.Assign(std::string("TC_0__1__34_2"),2,0,info);
  m_kZS_0.Assign(std::string("ZS_0"),2,0,info);
  p_vegas = new Vegas(m_rannum,100,m_name);
}

C3_13::~C3_13()
{
  delete p_vegas;
}

void C3_13::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 0;
  width = 0.;
}

void C3_13::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,p_rans);
}
std::string C3_13::ChID()
{
  return std::string("CGND$I_3_4$MTH_34$TC_0__1__34_2$ZS_0$");
}
