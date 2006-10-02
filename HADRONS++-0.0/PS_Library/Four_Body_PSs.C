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
	ResonanceFlavour prop1, 
	const int _k,
	ResonanceFlavour prop2,
	const int _i,
	const int _j ) :
  Single_Channel(1,4,fl),
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

void TwoResonances::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  PRINT_INFO("not implemented");
  abort();
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

void TwoResonances::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  PRINT_INFO("not implemented");
  abort();
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TwoResonancesParallel::TwoResonancesParallel( const Flavour * fl,
                                              ResonanceFlavour prop1,
                                              const int i,
                                              const int j,
                                              ResonanceFlavour prop2,
                                              const int k,
                                              const int l ) :
Single_Channel(1,4,fl),
m_P(Vec4D(fl[0].Mass(),0.,0.,0.)),
m_prop1 (prop1), m_prop2 (prop2),
m_I (i), m_J (j), m_K (k), m_L (l)
{
  name = string("TwoResonancesParallel_")
    + prop1.Name() + string("_")
    + ToString(m_I) + ToString(m_J)
    + string("_") + prop2.Name() + string("_")
    + ToString(m_K)+ToString(m_L);

  p_fl = new Flavour[5];
  for (short int i=0;i<nin+nout;i++) {
    p_fl[i] = fl[i];
    ms[i] = sqr(fl[i].PSMass());
  }
  msg.Tracking()<<"Init TwoResonancesParallel("<<name<<") : "<<endl
    <<"     "<<fl[0]<<" -> "
    <<fl[m_I]<<" "<<fl[m_J]<<" "<<fl[m_K]<<" "<<fl[m_L]<<", "<<endl
    <<"     "<<ms[0]<<" -> "
    <<ms[m_I]<<" "<<ms[m_J]<<" "<<ms[m_K]<<" "<<ms[m_L]<<endl
    <<"  => "<<p_fl[0]<<" -> "<<m_prop1.Name()<<" "<<m_prop2.Name()<<endl
    <<"     "<<p_fl[0]<<" -> "<<p_fl[m_I]<<" "<<p_fl[m_J]<<" "<<p_fl[m_K]<<" "<<p_fl[m_L]<<endl;

  rannum = 8;
  rans = new double[rannum];
  Integration_Info *info;
  info = new Integration_Info();
  m_kI_IJ_KL.Assign(std::string("I_IJ_KL"),2,0,info);
  m_kI_I_J.Assign(std::string("I_I_J"),2,0,info);
  m_kI_K_L.Assign(std::string("I_K_L"),2,0,info);
  m_kZRIJ_155.Assign(std::string("ZRIJ_155"),2,0,info);
  p_vegas = new Vegas(rannum,100,name);
}

void TwoResonancesParallel::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran)
{
  PRINT_INFO("not implemented");
  abort();
}

void TwoResonancesParallel::GeneratePoint(ATOOLS::Vec4D * p,double * _ran)
{
  // 0 -> prop1[-> I J] prop2[-> K L]
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<rannum;i++) rans[i]=ran[i];
  Vec4D  pIJKL = p[0];
  double sIJKL = dabs(pIJKL.Abs2());
  double sIJ_min = sqr( sqrt(ms[m_I]) + sqrt(ms[m_J]) );
  double sKL_min = sqr( sqrt(ms[m_K]) + sqrt(ms[m_L]) );
  double sKL_max = sqr(sqrt(sIJKL)-sqrt(sIJ_min));
  Vec4D  pKL;
  double sKL;
  sKL = CE.MassivePropMomenta(m_prop2.Mass(),m_prop2.Width(),1,sKL_min,sKL_max,ran[0]);
  double sIJ_max = sqr(sqrt(sIJKL)-sqrt(sKL));
  Vec4D  pIJ;
  double sIJ;
  sIJ = CE.MassivePropMomenta(m_prop1.Mass(),m_prop1.Width(),1,sIJ_min,sIJ_max,ran[1]);
  CE.Isotropic2Momenta(pIJKL,sIJ,sKL,pIJ,pKL,ran[2],ran[3]);
  double sJ = ms[m_J];
  double sI = ms[m_I];
  CE.Isotropic2Momenta(pIJ,sI,sJ,p[m_I],p[m_J],ran[4],ran[5]);
  double sK = ms[m_K];
  double sL = ms[m_L];
  CE.Isotropic2Momenta(pKL,sK,sL,p[m_K],p[m_L],ran[6],ran[7]);
}

void TwoResonancesParallel::GenerateWeight(Vec4D* p,Cut_Data * cuts)
{
  PRINT_INFO("not implemented");
  abort();
}

void TwoResonancesParallel::GenerateWeight(ATOOLS::Vec4D * p)
{
  // 0 -> prop1[-> I J] prop2[-> K L]
  double wt = 1.;
  Vec4D  pIJKL = p[0];
  double sIJKL = dabs(pIJKL.Abs2());
  double sIJ_min = sqr( sqrt(ms[m_I]) + sqrt(ms[m_J]) );
  double sKL_min = sqr( sqrt(ms[m_K]) + sqrt(ms[m_L]) );
  double sKL_max = sqr(sqrt(sIJKL)-sqrt(sIJ_min));
  Vec4D  pKL = p[m_K]+p[m_L];
  double sKL = dabs(pKL.Abs2());
  wt *= CE.MassivePropWeight(m_prop2.Mass(),m_prop2.Width(),1,sKL_min,sKL_max,sKL,rans[0]);
  double sIJ_max = sqr(sqrt(sIJKL)-sqrt(sKL));
  Vec4D  pIJ = p[m_J]+p[m_I];
  double sIJ = dabs(pIJ.Abs2());
  wt *= CE.MassivePropWeight(m_prop1.Mass(),m_prop1.Width(),1,sIJ_min,sIJ_max,sIJ,rans[1]);
  if (m_kI_IJ_KL.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_IJ_KL<<CE.Isotropic2Weight(pIJ,pKL,m_kI_IJ_KL[0],m_kI_IJ_KL[1]);
  wt *= m_kI_IJ_KL.Weight();
  
  rans[2]= m_kI_IJ_KL[0];
  rans[3]= m_kI_IJ_KL[1];
  if (m_kI_I_J.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_I_J<<CE.Isotropic2Weight(p[m_I],p[m_J],m_kI_I_J[0],m_kI_I_J[1]);
  wt *= m_kI_I_J.Weight();
  
  rans[4]= m_kI_I_J[0];
  rans[5]= m_kI_I_J[1];
  if (m_kI_K_L.Weight()==ATOOLS::UNDEFINED_WEIGHT)
    m_kI_K_L<<CE.Isotropic2Weight(p[m_K],p[m_L],m_kI_K_L[0],m_kI_K_L[1]);
  wt *= m_kI_K_L.Weight();
  
  rans[6]= m_kI_K_L[0];
  rans[7]= m_kI_K_L[1];
  double vw = p_vegas->GenerateWeight(rans);
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*3.-4.);
  
  weight = wt;
}
