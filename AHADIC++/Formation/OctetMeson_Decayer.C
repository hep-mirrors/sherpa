#include "AHADIC++/Formation/OctetMeson_Decayer.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

OctetMeson_Decayer::OctetMeson_Decayer(list<Singlet *> * singlets,
				       std::list<Proto_Particle *> * hadrons) :
  p_singlets(singlets), p_hadrons(hadrons),
  m_offset(Flavour(kf_eta_c_1S_oct).Kfcode()-Flavour(kf_eta_c_1S).Kfcode()),
  m_kappa(2.), m_minE(0.1) {}
    
bool OctetMeson_Decayer::operator()() {
  for (list<Singlet *>::iterator sit=p_singlets->begin();
       sit!=p_singlets->end();sit++) {
    p_part1 = p_part2 = NULL;
    //msg_Out()<<METHOD<<":\n"<<(**sit)<<"\n";
    if (!FixSpectator(*sit)) continue;
    if (!FixKinematics()) return false;
    UpdateColouredObjectsAndAddHadron();
  }
  return true;
}

bool OctetMeson_Decayer::FixSpectator(Singlet * singlet) {
  for (list<Proto_Particle *>::iterator ppiter1=singlet->begin();
       ppiter1!=singlet->end();ppiter1++) {
    p_part1 = (*ppiter1);
    //msg_Out()<<"Test "<<p_part1->Flavour()<<" --> "<<p_part1->Flavour().IsOctetMeson()<<"\n";
    if (p_part1->Flavour().IsOctetMeson()) {
      list<Proto_Particle *>::iterator ppiter2 = ppiter1, ppiter3; ppiter2++;
      if (ppiter2==singlet->end())     ppiter2 = p_singlet->begin();
      ppiter3 = (ppiter1==singlet->begin()) ? p_singlet->end() : ppiter1;
      ppiter3--;
      p_part2 = ( ( (p_part1->Momentum()+(*ppiter2)->Momentum()).Abs2() >
		    (p_part1->Momentum()+(*ppiter3)->Momentum()).Abs2() ) ?
		  (*ppiter2) : (*ppiter3) );
      break;
    }
  }
  return (p_part2!=NULL);
}

bool OctetMeson_Decayer::FixKinematics() {
  //msg_Out()<<METHOD<<"("<<p_part1<<"|"<<p_part2<<"):\n"
  //	   <<(*p_part1)
  //	   <<(*p_part2)<<"\n";
  // p+/- = E (1, 0, 0, +-1)
  // p1:    a1  p+ + (1-b1) p- --> m1^2/Q^2 = (1-b1) a1  ==> a1(1-b1) = m1^2/Q^2 
  // p2: (1-a1) p+ +    b1  p- --> m2^2/Q^2 = (1-a1) b1  ==> b1       = m2^2/Q^2 1/(1-a1)
  // ==> 0  = a1 [(1-a1)-m2^2/Q^2] - m1^2/Q^2 (1-a1)
  // ==> 0  = a1 [a1-1+m2^2/Q^2] - m1^2/Q^2 (a1-1) = 0
  // ==> 0  = Q^2 a1^2 - a1(Q^2+m1^2-m2^2) + m1^2
  // ==> a1 = 1/(2 Q^2) { Q^2+m1^2-m2^2 + sqrt[(Q^2-m1^2-m2^2)^2+4m2^2m1^2]}
  Vec4D  mom1      = p_part1->Momentum(), mom2 = p_part2->Momentum(), mom = mom1+mom2; 
  double Q2        = mom.Abs2(), E = sqrt(Q2)/2.;
  double m1        = p_part1->Flavour().HadMass(), m12 = sqr(m1);
  double m2        = p_part2->Flavour().HadMass(), m22 = sqr(m2);
  if (Q2<m12+m22) return false;
  double beta1     = 1./(2.*Q2) * (Q2-m12+m22 + sqrt(sqr(Q2-m12-m22)+4.*m12*m22));
  double alpha1    = m12/Q2/(1.-beta1);
  Poincare intoCMS = Poincare(mom);
  intoCMS.Boost(mom1);
  Poincare ontoZ   = Poincare(mom1,E*s_AxisP);
  ontoZ.Rotate(mom1);
  double zmin      = m_minE/(mom1[0]-m1), zmax = 1.-m12/(Q2-m22);
  if (zmax<zmin) zmin = zmax * m_minE/E;
  double z         = pow(pow(zmax,1.-m_kappa)+(1.-ran->Get())*pow(zmin,1.-m_kappa),
			 1./(1.-m_kappa));
  // p11: (1-z) a2  p+ + (1-b2) p- --> m1^2/Q^2 = (1-b2) (1-z) a2  ==> a2 = m1^2/Q^2 1/[(1-b2)(1-z)]
  // p12:    z  a2  p+             --> mg^2 = 0
  // p2:  (1   -a2) p+ +    b2  p1 --> m2^2/Q^2 = (1-a2) b2        ==> b2 = m2^2/Q^2 1/(1-a2) 
  // ==> 0 = a2^2 -a2(1-m1^2/Q^2/(1-z)-m2^2/Q^2 + m1^2/Q^2/(1-z)
  // have to check: (Q^2-m1^2/(1-z)-m2^2)^2 + 4m1^2 m2^2/(1-z) >
  // replace with Q^2-m2^2 > m1^2/(1-z) ==> 1-z > m1^2/(Q^2-m2^2) ==> z < zmax = 1-m1^2/(Q^2-m2^2) 
  m12 = m12/(1.-z);
  double beta2  = 1./(2.*Q2) * (Q2-m12+m22 + sqrt(sqr(Q2-m12-m22)+4.*m12*m22));
  double alpha2 = m12/Q2/(1.-beta2);
  m_mom[0]      = (1.-z)*alpha2 *E*s_AxisP + (1.-beta2)*E*s_AxisM;
  m_mom[1]      =     z *alpha2 *E*s_AxisP;
  m_mom[2]      = (1.-   alpha2)*E*s_AxisP +     beta2 *E*s_AxisM;
  for (size_t i=0;i<3;i++) {
    ontoZ.Rotate(m_mom[i]);
    intoCMS.BoostBack(m_mom[i]);
  }
  //msg_Out()<<"Selectec z = "<<z<<" in ["<<zmin<<", "<<zmax<<"]  --> "<<alpha2<<" & "<<beta2<<"\n"
  //	   <<p_part1->Momentum()<<" + "<<p_part2->Momentum()<<"\n"
  //	   <<m_mom[0]<<" + " <<m_mom[1]<<" + " <<m_mom[2]<<"\n"
  //	   <<"Check: "<<mom<<" vs. "<<(m_mom[0]+m_mom[1]+m_mom[2])<<"\n";
  return true;
}

void OctetMeson_Decayer::UpdateColouredObjectsAndAddHadron() {
  int newkfc = p_part1->Flavour().Kfcode()-m_offset;
  Proto_Particle * meson = new Proto_Particle(Flavour(newkfc),m_mom[0]);
  p_hadrons->push_back(meson);
  p_part1->SetFlavour(Flavour(kf_gluon));
  p_part1->SetMomentum(m_mom[1]);
  p_part2->SetMomentum(m_mom[2]);
  //msg_Out()<<(*meson)<<"\n"<<(*p_part1)<<"\n"<<(*p_part2)<<"\n";
}
