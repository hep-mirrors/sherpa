#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {
  class C4_16 : public Single_Channel {
    Info_Key m_kI_2_4,m_kI_3_24,m_kI_5_234,m_kZS_120;
    Vegas* p_vegas;
  public:
    C4_16(int,int,Flavour*,Integration_Info * const);
    ~C4_16();
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

extern "C" Single_Channel * Getter_C4_16(int nin,int nout,Flavour* fl,Integration_Info * const info) {
  return new C4_16(nin,nout,fl,info);
}

void C4_16::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D p2345=p[0]+p[1];
  double s2345_max = p2345.Abs2();
  double s234_max = sqr(sqrt(s2345_max)-sqrt(ms[5]));
  double s3 = ms[3];
  double s35_min = cuts->Getscut(std::string("35"));
  double s24_max = sqr(sqrt(s2345_max)-sqrt(s35_min));
  double s2 = ms[2];
  double s4 = ms[4];
  double s24_min = cuts->Getscut(std::string("24"));
  Vec4D  p24;
  double s24 = CE.ThresholdMomenta(1.5,4.*sqrt(s24_min),s24_min,s24_max,ran[0]);
  double s234_min = cuts->Getscut(std::string("234"));
  s234_min = Max(s234_min,sqr(sqrt(s24)+sqrt(s3)));
  Flavour fl234 = Flavour((kf_code)(24));
  Vec4D  p234;
  double s234 = CE.MassivePropMomenta(fl234.Mass(),fl234.Width(),1,s234_min,s234_max,ran[1]);
  double s5 = ms[5];
  CE.Isotropic2Momenta(p2345,s5,s234,p[5],p234,ran[2],ran[3]);
  CE.Isotropic2Momenta(p234,s3,s24,p[3],p24,ran[4],ran[5]);
  CE.Isotropic2Momenta(p24,s2,s4,p[2],p[4],ran[6],ran[7]);
}

void C4_16::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D p2345=p[0]+p[1];
  double s2345_max = p2345.Abs2();
  double s234_max = sqr(sqrt(s2345_max)-sqrt(ms[5]));
  double s3 = ms[3];
  double s35_min = cuts->Getscut(std::string("35"));
  double s24_max = sqr(sqrt(s2345_max)-sqrt(s35_min));
  double s2 = ms[2];
  double s4 = ms[4];
  double s24_min = cuts->Getscut(std::string("24"));
  Vec4D  p24 = p[2]+p[4];
  double s24 = dabs(p24.Abs2());
  wt *= CE.ThresholdWeight(1.5,4.*sqrt(s24_min),s24_min,s24_max,s24,rans[0]);
  double s234_min = cuts->Getscut(std::string("234"));
  s234_min = Max(s234_min,sqr(sqrt(s24)+sqrt(s3)));
  Flavour fl234 = Flavour((kf_code)(24));
  Vec4D  p234 = p[2]+p[3]+p[4];
  double s234 = dabs(p234.Abs2());
  wt *= CE.MassivePropWeight(fl234.Mass(),fl234.Width(),1,s234_min,s234_max,s234,rans[1]);
  double s5 = ms[5];
  if (m_kI_5_234.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_5_234<<CE.Isotropic2Weight(p[5],p234,m_kI_5_234[0],m_kI_5_234[1]);
  wt *= m_kI_5_234.Weight();

  rans[2]= m_kI_5_234[0];
  rans[3]= m_kI_5_234[1];
  if (m_kI_3_24.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_3_24<<CE.Isotropic2Weight(p[3],p24,m_kI_3_24[0],m_kI_3_24[1]);
  wt *= m_kI_3_24.Weight();

  rans[4]= m_kI_3_24[0];
  rans[5]= m_kI_3_24[1];
  if (m_kI_2_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_2_4<<CE.Isotropic2Weight(p[2],p[4],m_kI_2_4[0],m_kI_2_4[1]);
  wt *= m_kI_2_4.Weight();

  rans[6]= m_kI_2_4[0];
  rans[7]= m_kI_2_4[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*3.-4.);

  weight = wt;
}

C4_16::C4_16(int nin,int nout,Flavour* fl,Integration_Info * const info)
       : Single_Channel(nin,nout,fl)
{
  name = std::string("C4_16");
  rannum = 8;
  rans  = new double[rannum];
  m_kI_2_4.Assign(std::string("I_2_4"),2,0,info);
  m_kI_3_24.Assign(std::string("I_3_24"),2,0,info);
  m_kI_5_234.Assign(std::string("I_5_234"),2,0,info);
  m_kZS_120.Assign(std::string("ZS_120"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

C4_16::~C4_16()
{
  delete p_vegas;
}

void C4_16::ISRInfo(int & type,double & mass,double & width)
{
  type  = 2;
  mass  = 120.6285;
  width = 0.;
}

void C4_16::AddPoint(double Value)
{
  Single_Channel::AddPoint(Value);
  p_vegas->AddPoint(Value,rans);
}
std::string C4_16::ChID()
{
  return std::string("CG2$I_2_4$I_3_24$I_5_234$MP24_234$MTH_24$ZS_120$");
}
