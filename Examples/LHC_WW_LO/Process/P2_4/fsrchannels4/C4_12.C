#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C4_12 : public Single_Channel {
    double m_amct,m_alpha,m_ctmax,m_ctmin;
    Info_Key m_kI_2_5,m_kI_3_4,m_kTC_0__1__25_34,m_kZS_120;
    Vegas* p_vegas;
  public:
    C4_12(int,int,Flavour*,Integration_Info * const);
    ~C4_12();
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

extern "C" Single_Channel * Getter_C4_12(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C4_12(nin,nout,fl,info);
}

void C4_12::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p2345=p[0]+p[1];
  double s2345_max = p2345.Abs2();
  double s25_min = cuts->Getscut(std::string("25"));
  double s34_max = sqr(sqrt(s2345_max)-sqrt(s25_min));
  double s3 = ms[3];
  double s4 = ms[4];
  double s34_min = cuts->Getscut(std::string("34"));
  Flavour fl34 = Flavour((kf_code)(24));
  Vec4D  p34;
  double s34 = CE.MassivePropMomenta(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,ran[0]);
  s34_min = s34;
  double s25_max = sqr(sqrt(s2345_max)-sqrt(s34_min));
  double s2 = ms[2];
  double s5 = ms[5];
  Vec4D  p25;
  double s25 = CE.ThresholdMomenta(1.5,4.*sqrt(s25_min),s25_min,s25_max,ran[1]);
  CE.TChannelMomenta(p[0],p[1],p25,p34,s25,s34,0.,m_alpha,1.,-1.,m_amct,0,ran[2],ran[3]);
  CE.Isotropic2Momenta(p25,s2,s5,p[2],p[5],ran[4],ran[5]);
  CE.Isotropic2Momenta(p34,s3,s4,p[3],p[4],ran[6],ran[7]);
}

void C4_12::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p2345=p[0]+p[1];
  double s2345_max = p2345.Abs2();
  double s25_min = cuts->Getscut(std::string("25"));
  double s34_max = sqr(sqrt(s2345_max)-sqrt(s25_min));
  double s3 = ms[3];
  double s4 = ms[4];
  double s34_min = cuts->Getscut(std::string("34"));
  Flavour fl34 = Flavour((kf_code)(24));
  Vec4D  p34 = p[3]+p[4];
  double s34 = dabs(p34.Abs2());
  wt *= CE.MassivePropWeight(fl34.Mass(),fl34.Width(),1,s34_min,s34_max,s34,rans[0]);
  s34_min = s34;
  double s25_max = sqr(sqrt(s2345_max)-sqrt(s34_min));
  double s2 = ms[2];
  double s5 = ms[5];
  Vec4D  p25 = p[2]+p[5];
  double s25 = dabs(p25.Abs2());
  wt *= CE.ThresholdWeight(1.5,4.*sqrt(s25_min),s25_min,s25_max,s25,rans[1]);
  if (m_kTC_0__1__25_34.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kTC_0__1__25_34<<CE.TChannelWeight(p[0],p[1],p25,p34,0.,m_alpha,1.,-1.,m_amct,0,m_kTC_0__1__25_34[0],m_kTC_0__1__25_34[1]);
  wt *= m_kTC_0__1__25_34.Weight();

  rans[2]= m_kTC_0__1__25_34[0];
  rans[3]= m_kTC_0__1__25_34[1];
  if (m_kI_2_5.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_2_5<<CE.Isotropic2Weight(p[2],p[5],m_kI_2_5[0],m_kI_2_5[1]);
  wt *= m_kI_2_5.Weight();

  rans[4]= m_kI_2_5[0];
  rans[5]= m_kI_2_5[1];
  if (m_kI_3_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_3_4<<CE.Isotropic2Weight(p[3],p[4],m_kI_3_4[0],m_kI_3_4[1]);
  wt *= m_kI_3_4.Weight();

  rans[6]= m_kI_3_4[0];
  rans[7]= m_kI_3_4[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*3.-4.);

  weight = wt;
}

C4_12::C4_12(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C4_12");
  rannum = 8;
  rans  = new double[rannum];
  m_amct  = 1.;
  m_alpha = .9;
  m_ctmax = 1.;
  m_ctmin = -1.;
  m_kI_2_5.Assign(std::string("I_2_5"),2,0,info);
  m_kI_3_4.Assign(std::string("I_3_4"),2,0,info);
  m_kTC_0__1__25_34.Assign(std::string("TC_0__1__25_34"),2,0,info);
  m_kZS_120.Assign(std::string("ZS_120"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

C4_12::~C4_12()
{
  delete p_vegas;
}

void C4_12::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 120.6285;
  width = 0.;
}

void C4_12::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C4_12::ChID()
{
  return std::string("CG2$I_2_5$I_3_4$MP24_34$MTH_25$TC_0__1__25_34$ZS_120$");
}
