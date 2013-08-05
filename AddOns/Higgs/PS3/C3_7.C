#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C3_7 : public Single_Channel {
    Info_Key m_kI_2_34,m_kI_3_4,m_kZS_0;
    Vegas* p_vegas;
  public:
    C3_7(int,int,Flavour*,Integration_Info * const);
    ~C3_7();
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

extern "C" Single_Channel * Getter_C3_7(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C3_7(nin,nout,fl,info);
}

void C3_7::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s34_max = sqr(sqrt(s234_max)-sqrt(ms[2]));
  double s4 = ms[4];
  double s3 = ms[3];
  double s34_min = cuts->Getscut(std::string("34"));
  Vec4D  p34;
  double s34 = CE.MasslessPropMomenta(.5,s34_min,s34_max,ran[0]);
  double s2 = ms[2];
  CE.Isotropic2Momenta(p234,s2,s34,p[2],p34,ran[1],ran[2]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[3],ran[4]);
}

void C3_7::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s34_max = sqr(sqrt(s234_max)-sqrt(ms[2]));
  double s4 = ms[4];
  double s3 = ms[3];
  double s34_min = cuts->Getscut(std::string("34"));
  Vec4D  p34 = p[3]+p[4];
  double s34 = dabs(p34.Abs2());
  wt *= CE.MasslessPropWeight(.5,s34_min,s34_max,s34,rans[0]);
  double s2 = ms[2];
  if (m_kI_2_34.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_2_34<<CE.Isotropic2Weight(p[2],p34,m_kI_2_34[0],m_kI_2_34[1]);
  wt *= m_kI_2_34.Weight();

  rans[1]= m_kI_2_34[0];
  rans[2]= m_kI_2_34[1];
  if (m_kI_3_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_3_4<<CE.Isotropic2Weight(p[3],p[4],m_kI_3_4[0],m_kI_3_4[1]);
  wt *= m_kI_3_4.Weight();

  rans[3]= m_kI_3_4[0];
  rans[4]= m_kI_3_4[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,3*3.-4.);

  weight = wt;
}

C3_7::C3_7(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C3_7");
  rannum = 5;
  rans  = new double[rannum];
  m_kI_2_34.Assign(std::string("I_2_34"),2,0,info);
  m_kI_3_4.Assign(std::string("I_3_4"),2,0,info);
  m_kZS_0.Assign(std::string("ZS_0"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

C3_7::~C3_7()
{
  delete p_vegas;
}

void C3_7::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 0;
  width = 0.;
}

void C3_7::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C3_7::ChID()
{
  return std::string("CGND$I_2_34$I_3_4$MTH_34$ZS_0$");
}
