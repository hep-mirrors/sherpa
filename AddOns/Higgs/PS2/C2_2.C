#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C2_2 : public Single_Channel {
    Info_Key m_kI_2_3,m_kZR25_125;
    Vegas* p_vegas;
    int m_onshell;
  public:
    C2_2(int,int,Flavour*,Integration_Info * const);
    ~C2_2();
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

extern "C" Single_Channel * Getter_C2_2(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C2_2(nin,nout,fl,info);
}

void C2_2::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p23=p[0]+p[1];
  double s23_max = p23.Abs2();
  double s3 = ms[3];
  double s2 = ms[2];
  CE.Isotropic2Momenta(p23,s2,s3,p[2],p[3],ran[0],ran[1]);
}

void C2_2::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p23=p[0]+p[1];
  double s23_max = p23.Abs2();
  double s3 = ms[3];
  double s2 = ms[2];
  if (m_kI_2_3.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_2_3<<CE.Isotropic2Weight(p[2],p[3],m_kI_2_3[0],m_kI_2_3[1]);
  wt *= m_kI_2_3.Weight();

  rans[0]= m_kI_2_3[0];
  rans[1]= m_kI_2_3[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,2*3.-4.);

  weight = wt;
  if (m_onshell) {
    weight/=M_PI/p23.Abs2();
    weight*=Flavour((kf_code)(25)).Mass()*Flavour((kf_code)(25)).Width()*M_PI;
  }
}

C2_2::C2_2(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C2_2");
  rannum = 2;
  rans  = new double[rannum];
  m_kI_2_3.Assign(std::string("I_2_3"),2,0,info);
  m_kZR25_125.Assign(std::string("ZR25_125"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
  Data_Reader read(" ",";","#","=");
  m_onshell=read.GetValue<int>("HIGGS_ON_SHELL",0);
}

C2_2::~C2_2()
{
  delete p_vegas;
}

void C2_2::ISRInfo(int & type,double & mass,double & width)
{
  type  = 1;
  if (m_onshell) type=-1;
  mass  = Flavour((kf_code)(25)).Mass();
  width = Flavour((kf_code)(25)).Width();
}

void C2_2::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C2_2::ChID()
{
  return std::string("CG2$I_2_3$ZR25_125$");
}
