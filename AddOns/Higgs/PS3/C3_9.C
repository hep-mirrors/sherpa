#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C3_9 : public Single_Channel {
    Info_Key m_kI_2_4,m_kI_3_24,m_kZS_0;
    Vegas* p_vegas;
  public:
    C3_9(int,int,Flavour*,Integration_Info * const);
    ~C3_9();
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

extern "C" Single_Channel * Getter_C3_9(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C3_9(nin,nout,fl,info);
}

void C3_9::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s24_max = sqr(sqrt(s234_max)-sqrt(ms[3]));
  double s4 = ms[4];
  double s2 = ms[2];
  double s24_min = cuts->Getscut(std::string("24"));
  Vec4D  p24;
  double s24 = CE.MasslessPropMomenta(.5,s24_min,s24_max,ran[0]);
  double s3 = ms[3];
  CE.Isotropic2Momenta(p234,s3,s24,p[3],p24,ran[1],ran[2]);
  CE.Isotropic2Momenta(p24,s2,s4,p[2],p[4],ran[3],ran[4]);
}

void C3_9::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s24_max = sqr(sqrt(s234_max)-sqrt(ms[3]));
  double s24_min = cuts->Getscut(std::string("24"));
  Vec4D  p24 = p[2]+p[4];
  double s24 = dabs(p24.Abs2());
  wt *= CE.MasslessPropWeight(.5,s24_min,s24_max,s24,rans[0]);
  if (m_kI_3_24.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_3_24<<CE.Isotropic2Weight(p[3],p24,m_kI_3_24[0],m_kI_3_24[1]);
  wt *= m_kI_3_24.Weight();

  rans[1]= m_kI_3_24[0];
  rans[2]= m_kI_3_24[1];
  if (m_kI_2_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_2_4<<CE.Isotropic2Weight(p[2],p[4],m_kI_2_4[0],m_kI_2_4[1]);
  wt *= m_kI_2_4.Weight();

  rans[3]= m_kI_2_4[0];
  rans[4]= m_kI_2_4[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,3*3.-4.);

  weight = wt;
}

C3_9::C3_9(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C3_9");
  rannum = 5;
  rans  = new double[rannum];
  m_kI_2_4.Assign(std::string("I_2_4"),2,0,info);
  m_kI_3_24.Assign(std::string("I_3_24"),2,0,info);
  m_kZS_0.Assign(std::string("ZS_0"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

C3_9::~C3_9()
{
  delete p_vegas;
}

void C3_9::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 0;
  width = 0.;
}

void C3_9::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C3_9::ChID()
{
  return std::string("CGND$I_2_4$I_3_24$MTH_24$ZS_0$");
}
