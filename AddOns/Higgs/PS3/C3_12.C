#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C3_12 : public Single_Channel {
    double m_alpha,m_ctmax,m_ctmin;
    Info_Key m_kI_2_3,m_kTC_0__1__23_4,m_kZS_187;
    Vegas* p_vegas;
    bool m_onshell;
  public:
    C3_12(int,int,Flavour*,Integration_Info * const);
    ~C3_12();
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

extern "C" Single_Channel * Getter_C3_12(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C3_12(nin,nout,fl,info);
}

void C3_12::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<m_rannum;i++) p_rans[i]=ran[i];
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s4 = p_ms[4];
  double s23_max = sqr(sqrt(s234_max)-sqrt(p_ms[4]));
  double s3 = p_ms[3];
  double s2 = p_ms[2];
  double s23_min = cuts->Getscut((1<<2)|(1<<3));
  Flavour fl23 = Flavour((kf_code)(25));
  Vec4D  p23;
  double s23 = CE.MassivePropMomenta(fl23.Mass(),fl23.Width(),s23_min,s23_max,ran[0]);
  if (m_onshell) s23=sqr(fl23.Mass());
  m_ctmax = cuts->cosmax[1][4];
  m_ctmin = cuts->cosmin[1][4];
  CE.TChannelMomenta(p[0],p[1],p23,p[4],s23,s4,0.,m_alpha,m_ctmax,m_ctmin,ran[1-(size_t)m_onshell],ran[2-(size_t)m_onshell]);
  CE.Isotropic2Momenta(p23,s2,s3,p[2],p[3],ran[3-(size_t)m_onshell],ran[4-(size_t)m_onshell]);
}

void C3_12::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p234=p[0]+p[1];
  double s234_max = p234.Abs2();
  double s23_max = sqr(sqrt(s234_max)-sqrt(p_ms[4]));
  double s23_min = cuts->Getscut((1<<2)|(1<<3));
  Flavour fl23 = Flavour((kf_code)(25));
  Vec4D  p23 = p[2]+p[3];
  double s23 = dabs(p23.Abs2());
  if (m_onshell) wt /= (fl23.Mass()*fl23.Width())*M_PI;
  else
  wt *= CE.MassivePropWeight(fl23.Mass(),fl23.Width(),s23_min,s23_max,s23,p_rans[0]);
  m_ctmax = cuts->cosmax[1][4];
  m_ctmin = cuts->cosmin[1][4];
  if (m_kTC_0__1__23_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kTC_0__1__23_4<<CE.TChannelWeight(p[0],p[1],p23,p[4],0.,m_alpha,m_ctmax,m_ctmin,m_kTC_0__1__23_4[0],m_kTC_0__1__23_4[1]);
  wt *= m_kTC_0__1__23_4.Weight();

  p_rans[1-(size_t)m_onshell]= m_kTC_0__1__23_4[0];
  p_rans[2-(size_t)m_onshell]= m_kTC_0__1__23_4[1];
  if (m_kI_2_3.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_2_3<<CE.Isotropic2Weight(p[2],p[3],m_kI_2_3[0],m_kI_2_3[1]);
  wt *= m_kI_2_3.Weight();

  p_rans[3-(size_t)m_onshell]= m_kI_2_3[0];
  p_rans[4-(size_t)m_onshell]= m_kI_2_3[1];
  double vw = p_vegas->GenerateWeight(p_rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,3*3.-4.);

  m_weight = wt;
}

C3_12::C3_12(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  Settings& s = Settings::GetMainSettings();
  m_onshell = s["HIGGS_ON_SHELL"].Get<bool>();
  m_name = std::string("C3_12");
  m_rannum = 5-(size_t)m_onshell;
  p_rans  = new double[m_rannum];
  m_alpha = s["SCHANNEL_ALPHA"].Get<double>();
  m_ctmax = 1.;
  m_ctmin = -1.;
  m_kI_2_3.Assign(std::string("I_2_3"),2,0,info);
  m_kTC_0__1__23_4.Assign(std::string("TC_0__1__23_4"),2,0,info);
  m_kZS_187.Assign(std::string("ZS_187"),2,0,info);
  p_vegas = new Vegas(m_rannum,100,m_name);
}

C3_12::~C3_12()
{
  delete p_vegas;
}

void C3_12::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 187.5;
  width = 0.;
}

void C3_12::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,p_rans);
}
std::string C3_12::ChID()
{
  return std::string("CGND$I_2_3$MP25_23$TC_0__1__23_4$ZS_187$");
}
