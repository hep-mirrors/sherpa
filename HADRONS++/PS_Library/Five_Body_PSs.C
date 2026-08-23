#include "HADRONS++/PS_Library/Five_Body_PSs.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"


using namespace HADRONS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// ThreeResonances: direct extension of Four_Body_PSs.C's TwoResonances
// by one more sequential resonance. Every formula below mirrors
// TwoResonances' own structure one-for-one, just with an extra
// (mass,index) pair (m_l, m_prop2 moved down one level, m_prop3 added
// as the new innermost resonance) - see that file for the general
// design rationale (Breit-Wigner-peaked mass sampling at each level,
// nested 2-body isotropic splits, Jacobian-matched weight in
// GenerateWeight()).
//
///////////////////////////////////////////////////////////////////////////

ThreeResonances::ThreeResonances(const Flavour * fl,
				  SimpleResonanceFlavour prop1,const int _k,
				  SimpleResonanceFlavour prop2,const int _l,
				  SimpleResonanceFlavour prop3,const int _i,
				  const int _j ) :
  Single_Channel(1,5,fl),
  m_P(Vec4D(fl[0].HadMass(),0.,0.,0.)),
  m_i (_i), m_j (_j), m_k (_k), m_l (_l),
  m_prop1 (prop1), m_prop2 (prop2), m_prop3 (prop3)
{
  m_name = string("ThreeResonances_")
    + prop1.Name() + string("_")
    + ToString(m_k)
    + string("_") + prop2.Name() + string("_")
    + ToString(m_l)
    + string("_") + prop3.Name() + string("_")
    + ToString(m_i)+ToString(m_j);
  // generate channel name
  p_fl = new Flavour[6];
  for (short int i=0;i<m_nin+m_nout;i++) {
    p_fl[i] = fl[i];
    p_ms[i] = sqr(fl[i].HadMass());
  }
  // set masses^2
  for (int i=1;i<6;i++) {
    if (m_i!=i && m_j!=i && m_k!=i && m_l!=i) { m_dir=i; break; }
  }				// get the one with no resonance
  msg_Tracking()<<"Init ThreeResonances("<<m_name<<") : "<<endl
		<<"     "<<fl[0]<<" -> "
		<<fl[m_dir]<<" "<<fl[m_k]<<" "<<fl[m_l]<<" "<<fl[m_i]<<" "<<fl[m_j]
		<<", "<<endl
		<<"     "<<p_ms[0]<<" -> "
		<<p_ms[m_dir]<<" "<<p_ms[m_k]<<" "<<p_ms[m_l]<<" "<<p_ms[m_i]<<" "<<p_ms[m_j]<<endl
		<<"  => "<<p_fl[0]<<" -> "<<p_fl[m_dir]<<" "<<m_prop1.Name()
		<<endl
		<<"     "<<p_fl[0]<<" -> "<<p_fl[m_dir]<<" "<<p_fl[m_k]<<" "
		<<m_prop2.Name()<<endl
		<<"     "<<p_fl[0]<<" -> "<<p_fl[m_dir]<<" "<<p_fl[m_k]<<" "
		<<p_fl[m_l]<<" "<<m_prop3.Name()<<endl
		<<"     "<<p_fl[0]<<" -> "<<p_fl[m_dir]<<" "<<p_fl[m_k]<<" "
		<<p_fl[m_l]<<" "<<p_fl[m_i]<<" "<<p_fl[m_j]<<endl;
  msg_Debugging()
    <<"  with res1 @ "<<m_prop1.Mass()<<" ("<<m_prop1.Width()<<")"<<endl
    <<"      res2 @ "<<m_prop2.Mass()<<" ("<<m_prop2.Width()<<")"<<endl
    <<"      res3 @ "<<m_prop3.Mass()<<" ("<<m_prop3.Width()<<")"<<endl;

  m_rannum = 11;
  p_rans = new double[m_rannum];
  p_vegas = new Vegas(m_rannum,100,m_name);
  p_info = new Integration_Info();
  m_kI_1234_5.Assign(std::string("I_1234_5"),2,0,p_info);
  m_kI_123_4.Assign(std::string("I_123_4"),2,0,p_info);
  m_kI_12_3.Assign(std::string("I_12_3"),2,0,p_info);
  m_kI_1_2.Assign(std::string("I_1_2"),2,0,p_info);
}

ThreeResonances::~ThreeResonances() {
  if (p_fl!=NULL) delete [] p_fl; p_fl = NULL;
  if (p_vegas) delete p_vegas; p_vegas=NULL;
  if (p_info) delete p_info; p_info=NULL;
}

void ThreeResonances::GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data * cuts,
				     double * _ran)
{
  double *ran = p_vegas->GeneratePoint(_ran);
  for(int i=0;i<m_rannum;i++) p_rans[i]=ran[i];
  Vec4D  p12345 = p[0];
  // kinematic variables - innermost pair (i,j) outward to (k) and the
  // outer spectator (dir), building up minimal invariant masses level
  // by level exactly as TwoResonances does, just with one more level
  // (the new m_l/m_prop2 pairing sits between the old k/res2 and i/j).
  double s1_min = p_ms[m_i];
  double s2_min = p_ms[m_j];
  double s3_min = p_ms[m_l];
  double s4_min = p_ms[m_k];
  double s12_min  = sqr( sqrt(s1_min) + sqrt(s2_min) );
  double s123_min = sqr( sqrt(s12_min) + sqrt(s3_min) );
  double s1234_min = sqr( sqrt(s123_min) + sqrt(s4_min) );
  double s12345 = dabs(p12345.Abs2());
  double s1 = p_ms[m_i];
  double s2 = p_ms[m_j];
  double s3 = p_ms[m_l];
  double s4 = p_ms[m_k];
  double s5 = p_ms[m_dir];
  double s1234_max = sqr(sqrt(s12345)-sqrt(s5));
  Vec4D  p1234;
  double s1234;
  s1234 = CE.MassivePropMomenta(m_prop1.Mass(),m_prop1.Width(),
				 s1234_min,s1234_max,ran[0]);
  CE.Isotropic2Momenta(p12345,s1234,s5,p1234,p[m_dir],ran[1],ran[2]);
  double s123_max = sqr(sqrt(s1234)-sqrt(s4));
  Vec4D  p123;
  double s123;
  s123 = CE.MassivePropMomenta(m_prop2.Mass(),m_prop2.Width(),
				s123_min,s123_max,ran[3]);
  CE.Isotropic2Momenta(p1234,s123,s4,p123,p[m_k],ran[4],ran[5]);
  double s12_max = sqr(sqrt(s123)-sqrt(s3));
  Vec4D  p12;
  double s12;
  s12 = CE.MassivePropMomenta(m_prop3.Mass(),m_prop3.Width(),
			       s12_min,s12_max,ran[6]);
  CE.Isotropic2Momenta(p123,s12,s3,p12,p[m_l],ran[7],ran[8]);
  CE.Isotropic2Momenta(p12,s1,s2,p[m_i],p[m_j],ran[9],ran[10]);
}

void ThreeResonances::GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data * cuts)
{
  double wt = 1.;
  Vec4D  p12345 = p[0];
  double s1_min = p_ms[m_i];
  double s2_min = p_ms[m_j];
  double s3_min = p_ms[m_l];
  double s4_min = p_ms[m_k];
  double s12_min  = sqr( sqrt(s1_min) + sqrt(s2_min) );
  double s123_min = sqr( sqrt(s12_min) + sqrt(s3_min) );
  double s1234_min = sqr( sqrt(s123_min) + sqrt(s4_min) );
  double s12345 = dabs(p12345.Abs2());
  double s4 = p_ms[m_k];
  double s5 = p_ms[m_dir];
  double s1234_max = sqr(sqrt(s12345)-sqrt(s5));
  Vec4D  p1234 = p[m_i]+p[m_j]+p[m_l]+p[m_k];
  double s1234 = dabs(p1234.Abs2());
  wt *= CE.MassivePropWeight(m_prop1.Mass(),m_prop1.Width(),
			      s1234_min,s1234_max,s1234,p_rans[0]);
  m_kI_1234_5<<CE.Isotropic2Weight(p[0]-p[m_dir],p[m_dir],m_kI_1234_5[0],m_kI_1234_5[1]);
  wt *= m_kI_1234_5.Weight();

  p_rans[1]= m_kI_1234_5[0];
  p_rans[2]= m_kI_1234_5[1];
  double s3 = p_ms[m_l];
  double s123_max = sqr(sqrt(s1234)-sqrt(s4));
  Vec4D  p123 = p[m_i]+p[m_j]+p[m_l];
  double s123 = dabs(p123.Abs2());
  wt *= CE.MassivePropWeight(m_prop2.Mass(),m_prop2.Width(),
			      s123_min,s123_max,s123,p_rans[3]);
  m_kI_123_4<<CE.Isotropic2Weight(p123,p[m_k],m_kI_123_4[0],m_kI_123_4[1]);
  wt *= m_kI_123_4.Weight();

  p_rans[4]= m_kI_123_4[0];
  p_rans[5]= m_kI_123_4[1];
  double s12_max = sqr(sqrt(s123)-sqrt(s3));
  Vec4D  p12 = p[m_i]+p[m_j];
  double s12 = dabs(p12.Abs2());
  wt *= CE.MassivePropWeight(m_prop3.Mass(),m_prop3.Width(),
			      s12_min,s12_max,s12,p_rans[6]);
  m_kI_12_3<<CE.Isotropic2Weight(p12,p[m_l],m_kI_12_3[0],m_kI_12_3[1]);
  wt *= m_kI_12_3.Weight();

  p_rans[7]= m_kI_12_3[0];
  p_rans[8]= m_kI_12_3[1];
  m_kI_1_2<<CE.Isotropic2Weight(p[m_i],p[m_j],m_kI_1_2[0],m_kI_1_2[1]);
  wt *= m_kI_1_2.Weight();

  p_rans[9] = m_kI_1_2[0];
  p_rans[10]= m_kI_1_2[1];
  double vw = p_vegas->GenerateWeight(p_rans);
  // Phase-space volume normalisation: TwoResonances uses
  // pow(2 pi, 4*3-4) for nout=4 (i.e. 4*(nout-1)-4). Following the
  // same pattern for nout=5: 4*(5-1)-4 = 12.
  if (wt!=0.) wt = vw/wt/pow(2.*M_PI,4*4.-4.);

  m_weight = wt;
}
