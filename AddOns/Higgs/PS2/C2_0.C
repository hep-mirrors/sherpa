#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C2_0 : public Single_Channel {
    double m_amct,m_alpha,m_ctmax,m_ctmin;
    Info_Key m_kTC_0__1__2_3,m_kZS_0;
    Vegas* p_vegas;
  public:
    C2_0(int,int,Flavour*,Integration_Info * const);
    ~C2_0();
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

extern "C" Single_Channel * Getter_C2_0(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C2_0(nin,nout,fl,info);
}

void C2_0::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p23=p[0]+p[1];
  double s3 = ms[3];
  double s2 = ms[2];
  m_ctmax = Min(cuts->cosmax[0][2],cuts->cosmax[1][3]);
  CE.TChannelMomenta(p[0],p[1],p[2],p[3],s2,s3,0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,ran[0],ran[1]);
}

void C2_0::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p23=p[0]+p[1];
  m_ctmax = Min(cuts->cosmax[0][2],cuts->cosmax[1][3]);
  if (m_kTC_0__1__2_3.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kTC_0__1__2_3<<CE.TChannelWeight(p[0],p[1],p[2],p[3],0.,m_alpha,m_ctmax,m_ctmin,m_amct,0,m_kTC_0__1__2_3[0],m_kTC_0__1__2_3[1]);
  wt *= m_kTC_0__1__2_3.Weight();

  rans[0]= m_kTC_0__1__2_3[0];
  rans[1]= m_kTC_0__1__2_3[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,2*3.-4.);

  weight = wt;
}

C2_0::C2_0(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C2_0");
  rannum = 2;
  rans  = new double[rannum];
  m_amct  = 1.;
  m_alpha = .9;
  m_ctmax = 1.;
  m_ctmin = -1.;
  m_kTC_0__1__2_3.Assign(std::string("TC_0__1__2_3"),2,0,info);
  m_kZS_0.Assign(std::string("ZS_0"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

C2_0::~C2_0()
{
  delete p_vegas;
}

void C2_0::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 0;
  width = 0.;
}

void C2_0::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C2_0::ChID()
{
  return std::string("CG2$TC_0__1__2_3$ZS_0$");
}
