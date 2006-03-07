#include "Four_Body_PSs.H"
#include "Channel_Elements.H"
#include "Channel_Basics.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Random.H"


using namespace HADRONS; 
using namespace PHASIC; 
using namespace ATOOLS; 
using namespace std; 

TwoResonances::TwoResonances(
	const Flavour * fl,
	SimpleResonanceFlavour prop1, 
	const int _k,
	SimpleResonanceFlavour prop2,
	const int _i,
	const int _j )
: Single_Channel(1,4,fl), 
  m_P(Vec4D(fl[0].Mass(),0.,0.,0.)), 
  m_prop1 (prop1), m_prop2 (prop2),
  m_i (_i), m_j (_j), m_k (_k)
{
  name = string("TwoResonances_")
	+ prop1.Name() + string("_")
	+ ToString(m_k)
	+ string("_") + prop2.Name() + string("_")
	+ ToString(m_i)+ToString(m_j);
											// generate channel name
  p_fl = new Flavour[5];
  for (short int i=0;i<nin+nout;i++) {
	p_fl[i] = fl[i];
	ms[i] = sqr(fl[i].PSMass());
  }
											// set masses^2
  for (int i=1;i<5;i++) {
    if (m_i!=i && m_j!=i && m_k!=i) { m_dir=i; break; }
  }				// get the one with no resonance
  msg.Tracking()<<"Init TwoResonances("<<name<<") : "<<endl
	   <<"     "<<fl[0]<<" -> "
	            <<fl[m_dir]<<" "<<fl[m_k]<<" "<<fl[m_i]<<" "<<fl[m_j]<<", "<<endl
	   <<"     "<<ms[0]<<" -> "
	            <<ms[m_dir]<<" "<<ms[m_k]<<" "<<ms[m_i]<<" "<<ms[m_j]<<endl
	   <<"  => "<<p_fl[0]<<" -> "<<p_fl[m_dir]<<" "<<m_prop1.Name()<<endl
	   <<"     "<<p_fl[0]<<" -> "<<p_fl[m_dir]<<" "<<p_fl[m_k]<<" "<<m_prop2.Name()<<endl
	   <<"     "<<p_fl[0]<<" -> "<<p_fl[m_dir]<<" "<<p_fl[m_k]<<" "<<p_fl[m_i]<<" "<<p_fl[m_j]<<endl;
  msg.Debugging()
       <<"  with axial @ "<<m_prop1.Mass()<<" ("<<m_prop1.Width()<<")"<<endl
	   <<"      vector @ "<<m_prop2.Mass()<<" ("<<m_prop2.Width()<<")"<<endl;

  rannum = 8;
  rans = new double[rannum];
  p_vegas = new Vegas(rannum,100,name);
  Integration_Info *info;
  info = new Integration_Info();
  m_kI_123_4.Assign(std::string("I_123_4"),2,0,info);
  m_kI_12_3.Assign(std::string("I_12_3"),2,0,info);
  m_kI_1_2.Assign(std::string("I_1_2"),2,0,info);
}

void TwoResonances::GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * cuts,double * _ran)
{
  /*
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D  p1234 = p[0];
  double s1234 = dabs(p1234.Abs2());
  double s123_min = cuts->Getscut(std::string("123"));
  double s4 = ms[4];
  double s4_min = cuts->Getscut(std::string("4"));
  double s123_max = sqr(sqrt(s1234)-sqrt(s4_min));
  Flavour fl123 = Flavour(kf::code(20113));
  Vec4D  p123;
  double s123;
  s123 = CE.MassivePropMomenta(fl123.PSMass(),fl123.Width(),1,s123_min,s123_max,ran[0]);
  CE.Isotropic2Momenta(p1234,s123,s4,p123,p[4],ran[1],ran[2]);
  double s12_min = cuts->Getscut(std::string("12"));
  double s3 = ms[3];
  double s3_min = cuts->Getscut(std::string("3"));
  double s12_max = sqr(sqrt(s123)-sqrt(s3_min));
  Flavour fl12 = Flavour(kf::code(113));
  Vec4D  p12;
  double s12;
  s12 = CE.MassivePropMomenta(fl12.PSMass(),fl12.Width(),1,s12_min,s12_max,ran[3]);
  CE.Isotropic2Momenta(p123,s12,s3,p12,p[3],ran[4],ran[5]);
  double s1 = ms[1];
  double s2 = ms[2];
  CE.Isotropic2Momenta(p12,s1,s2,p[1],p[2],ran[6],ran[7]);
  */
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D  p1234 = p[0];
  // kinematic variables
  double s1234_min = ms[0];
  double s1_min = ms[m_i];
  double s2_min = ms[m_j];
  double s3_min = ms[m_k];
  double s4_min = ms[m_dir];
  double s12_min = sqr( sqrt(s1_min) + sqrt(s2_min) );
  double s123_min = sqr( sqrt(s12_min) + sqrt(s3_min) );
  double s1234 = dabs(p1234.Abs2());
  double s1 = ms[m_i];
  double s2 = ms[m_j];
  double s3 = ms[m_k];
  double s4 = ms[m_dir];
  double s123_max = sqr(sqrt(s1234)-sqrt(s4));
  Vec4D  p123;
  double s123;
  s123 = CE.MassivePropMomenta(m_prop1.Mass(),m_prop1.Width(),1,s123_min,s123_max,ran[0]);
  CE.Isotropic2Momenta(p1234,s123,s4,p123,p[m_dir],ran[1],ran[2]);
  double s12_max = sqr(sqrt(s123)-sqrt(s3));
  Vec4D  p12;
  double s12;
  s12 = CE.MassivePropMomenta(m_prop2.Mass(),m_prop2.Width(),1,s12_min,s12_max,ran[3]);
  CE.Isotropic2Momenta(p123,s12,s3,p12,p[m_k],ran[4],ran[5]);
  CE.Isotropic2Momenta(p12,s1,s2,p[m_i],p[m_j],ran[6],ran[7]);
}

void TwoResonances::GeneratePoint(ATOOLS::Vec4D * p,double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D  p1234 = p[0];
  // kinematic variables
  double s1234_min = ms[0];
  double s1_min = ms[m_i];
  double s2_min = ms[m_j];
  double s3_min = ms[m_k];
  double s4_min = ms[m_dir];
  double s12_min = sqr( sqrt(s1_min) + sqrt(s2_min) );
  double s123_min = sqr( sqrt(s12_min) + sqrt(s3_min) );
  double s1234 = dabs(p1234.Abs2());
  double s1 = ms[m_i];
  double s2 = ms[m_j];
  double s3 = ms[m_k];
  double s4 = ms[m_dir];
  double s123_max = sqr(sqrt(s1234)-sqrt(s4));
  Vec4D  p123;
  double s123;
  s123 = CE.MassivePropMomenta(m_prop1.Mass(),m_prop1.Width(),1,s123_min,s123_max,ran[0]);
  CE.Isotropic2Momenta(p1234,s123,s4,p123,p[m_dir],ran[1],ran[2]);
  double s12_max = sqr(sqrt(s123)-sqrt(s3));
  Vec4D  p12;
  double s12;
  s12 = CE.MassivePropMomenta(m_prop2.Mass(),m_prop2.Width(),1,s12_min,s12_max,ran[3]);
  CE.Isotropic2Momenta(p123,s12,s3,p12,p[m_k],ran[4],ran[5]);
  CE.Isotropic2Momenta(p12,s1,s2,p[m_i],p[m_j],ran[6],ran[7]);
}


void TwoResonances::GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * cuts)
{
  /*
  double wt = 1.;
  Vec4D  p1234 = p[0];
  double s1234 = dabs(p1234.Abs2());
  double s123_min = cuts->Getscut(std::string("123"));
  double s4_min = cuts->Getscut(std::string("4"));
  double s123_max = sqr(sqrt(s1234)-sqrt(s4_min));
  Flavour fl123 = Flavour(kf::code(20113));
  Vec4D  p123 = p[1]+p[2]+p[3];
  double s123 = dabs(p123.Abs2());
  wt *= CE.MassivePropWeight(fl123.PSMass(),fl123.Width(),1,s123_min,s123_max,s123,rans[0]);
  if (m_kI_123_4.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_123_4<<CE.Isotropic2Weight(p123,p[4],m_kI_123_4[0],m_kI_123_4[1]);
  wt *= m_kI_123_4.Weight();

  rans[1]= m_kI_123_4[0];
  rans[2]= m_kI_123_4[1];
  double s12_min = cuts->Getscut(std::string("12"));
  double s3_min = cuts->Getscut(std::string("3"));
  double s12_max = sqr(sqrt(s123)-sqrt(s3_min));
  Flavour fl12 = Flavour(kf::code(113));
  Vec4D  p12 = p[1]+p[2];
  double s12 = dabs(p12.Abs2());
  wt *= CE.MassivePropWeight(fl12.PSMass(),fl12.Width(),1,s12_min,s12_max,s12,rans[3]);
  if (m_kI_12_3.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_12_3<<CE.Isotropic2Weight(p12,p[3],m_kI_12_3[0],m_kI_12_3[1]);
  wt *= m_kI_12_3.Weight();

  rans[4]= m_kI_12_3[0];
  rans[5]= m_kI_12_3[1];
  if (m_kI_1_2.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_1_2<<CE.Isotropic2Weight(p[1],p[2],m_kI_1_2[0],m_kI_1_2[1]);
  wt *= m_kI_1_2.Weight();

  rans[6]= m_kI_1_2[0];
  rans[7]= m_kI_1_2[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*3.-4.);

  weight = wt;
  */
  double wt = 1.;
  Vec4D  p1234 = p[0];
  // kinematic variables
  double s1234_min = ms[0];
  double s1_min = ms[m_i];
  double s2_min = ms[m_j];
  double s3_min = ms[m_k];
  double s4_min = ms[m_dir];
  double s12_min = sqr( sqrt(s1_min) + sqrt(s2_min) );
  double s123_min = sqr( sqrt(s12_min) + sqrt(s3_min) );
  double s1234 = dabs(p1234.Abs2());
  double s1 = ms[m_i];
  double s2 = ms[m_j];
  double s3 = ms[m_k];
  double s4 = ms[m_dir];
  double s123_max = sqr(sqrt(s1234)-sqrt(s4));
  Vec4D  p123 = p[m_i]+p[m_j]+p[m_k];
  double s123 = dabs(p123.Abs2());
  wt *= CE.MassivePropWeight(m_prop1.Mass(),m_prop1.Width(),1,s123_min,s123_max,s123,rans[0]);
  m_kI_123_4<<CE.Isotropic2Weight(p123,p[m_dir],m_kI_123_4[0],m_kI_123_4[1]);
  wt *= m_kI_123_4.Weight();

  rans[1]= m_kI_123_4[0];
  rans[2]= m_kI_123_4[1];
  double s12_max = sqr(sqrt(s123)-sqrt(s3));
  Vec4D  p12 = p[m_i]+p[m_j];
  double s12 = dabs(p12.Abs2());
  wt *= CE.MassivePropWeight(m_prop2.Mass(),m_prop2.Width(),1,s12_min,s12_max,s12,rans[3]);
  m_kI_12_3<<CE.Isotropic2Weight(p12,p[m_k],m_kI_12_3[0],m_kI_12_3[1]);
  wt *= m_kI_12_3.Weight();
 
  rans[4]= m_kI_12_3[0];
  rans[5]= m_kI_12_3[1];
  m_kI_1_2<<CE.Isotropic2Weight(p[m_i],p[m_j],m_kI_1_2[0],m_kI_1_2[1]);
  wt *= m_kI_1_2.Weight();
 
  rans[6]= m_kI_1_2[0];
  rans[7]= m_kI_1_2[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*3.-4.);

  weight = wt;
}


void TwoResonances::GenerateWeight(ATOOLS::Vec4D * p)
{
  double wt = 1.;
  Vec4D  p1234 = p[0];
  // kinematic variables
  double s1234_min = ms[0];
  double s1_min = ms[m_i];
  double s2_min = ms[m_j];
  double s3_min = ms[m_k];
  double s4_min = ms[m_dir];
  double s12_min = sqr( sqrt(s1_min) + sqrt(s2_min) );
  double s123_min = sqr( sqrt(s12_min) + sqrt(s3_min) );
  double s1234 = dabs(p1234.Abs2());
  double s1 = ms[m_i];
  double s2 = ms[m_j];
  double s3 = ms[m_k];
  double s4 = ms[m_dir];
  double s123_max = sqr(sqrt(s1234)-sqrt(s4));
  Vec4D  p123 = p[m_i]+p[m_j]+p[m_k];
  double s123 = dabs(p123.Abs2());
  wt *= CE.MassivePropWeight(m_prop1.Mass(),m_prop1.Width(),1,s123_min,s123_max,s123,rans[0]);
  m_kI_123_4<<CE.Isotropic2Weight(p123,p[m_dir],m_kI_123_4[0],m_kI_123_4[1]);
  wt *= m_kI_123_4.Weight();

  rans[1]= m_kI_123_4[0];
  rans[2]= m_kI_123_4[1];
  double s12_max = sqr(sqrt(s123)-sqrt(s3));
  Vec4D  p12 = p[m_i]+p[m_j];
  double s12 = dabs(p12.Abs2());
  wt *= CE.MassivePropWeight(m_prop2.Mass(),m_prop2.Width(),1,s12_min,s12_max,s12,rans[3]);
  m_kI_12_3<<CE.Isotropic2Weight(p12,p[m_k],m_kI_12_3[0],m_kI_12_3[1]);
  wt *= m_kI_12_3.Weight();
 
  rans[4]= m_kI_12_3[0];
  rans[5]= m_kI_12_3[1];
  m_kI_1_2<<CE.Isotropic2Weight(p[m_i],p[m_j],m_kI_1_2[0],m_kI_1_2[1]);
  wt *= m_kI_1_2.Weight();
 
  rans[6]= m_kI_1_2[0];
  rans[7]= m_kI_1_2[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*3.-4.);

  weight = wt;
}
