#include "AHADIC++/Formation/OctetMeson_Decayer.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/LDME.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

OctetMeson_Decayer::OctetMeson_Decayer(list<Singlet *> *singlets,
                                       std::list<Proto_Particle *> *hadrons) :
  p_singlets(singlets), p_hadrons(hadrons),
  m_offset(Flavour(kf_1S0_c).Kfcode() -
	  Flavour(kf_eta_c_1S).Kfcode()),
  m_kappa(2.), m_minE(0.1) {}

bool OctetMeson_Decayer::operator()() {
  //msg_Out()<<METHOD<<" for "<<p_singlets->size()<<" singlets.\n";
  for (list<Singlet *>::iterator sit = p_singlets->begin();
       sit != p_singlets->end(); sit++) {
    p_part1 = p_part2 = NULL;
    list<Proto_Particle *>::iterator pit1;
    bool found;
    do {
      found = false;
      if (FindOctet(*sit, pit1) && FixSpectator(*sit, pit1)) {
        if (!FixKinematics())
          return false;
        UpdateColouredObjectsAndAddHadron();
        found = true;
      }
    } while (found);
    //msg_Out()<<(**sit)<<"\n";
  }
  return true;
}

bool OctetMeson_Decayer::FindOctet(Singlet *singlet,
                                   list<Proto_Particle *>::iterator &pit1) {
  for (pit1 = singlet->begin(); pit1 != singlet->end(); pit1++) {
    // msg_Out()<<"kf_code: "<<(*pit1)-> <<'\n';
    if ((*pit1)->Flavour().IsOctetMeson()) {
      p_part1 = (*pit1);
      kf_code octet_kfc = p_part1->Flavour().Kfcode();
      newkfc = GetQuarkoniumDecayChannel(octet_kfc);
      if(newkfc == octet_kfc) break; //decay did not happen
      return true;
    }
  }
  return false;
}

bool OctetMeson_Decayer::FixSpectator(Singlet *singlet,
                                      list<Proto_Particle *>::iterator &pit1) {
  list<Proto_Particle *>::iterator pit2 = pit1, pit3 = pit1;
  Proto_Particle *part2 = NULL, *part3 = NULL;
  if (pit2 != singlet->begin()) {
    pit2--;
    part2 = (*pit2);
  }
  if ((++pit3) != singlet->end()) {
    part3 = (*pit3);
  }
  if (part2 == NULL && part3 == NULL) {
    THROW(fatal_error, "this seems to be a one-particle singlet!");
  }
  else if (part2 == NULL)
    p_part2 = part3;
  else if (part3 == NULL)
    p_part2 = part2;
  else if (!part2->Flavour().IsOctetMeson() &&
	   part3->Flavour().IsOctetMeson())
    p_part2 = part2;
  else if (part2->Flavour().IsOctetMeson() &&
	   !part3->Flavour().IsOctetMeson())
    p_part2 = part3;
  else if (( (p_part1->Momentum() + part2->Momentum()).Abs2() -
	     p_part1->Momentum().Abs2()-part2->Momentum().Abs2() ) >
           ( (p_part1->Momentum() + part3->Momentum()).Abs2() -
	     p_part1->Momentum().Abs2()-part3->Momentum().Abs2() )) 
    p_part2 = part2;
  else
    p_part2 = part3;
  return true;
}

bool OctetMeson_Decayer::FixKinematics() {
  //msg_Out()<<"\n"<<METHOD<<"("<<p_part1<<"|"<<p_part2<<"|):\n"
  //	   <<(*p_part1)<<(*p_part2);
  //  p+/- = E (1, 0, 0, +-1)
  //  p1:    a1  p+ + (1-b1) p- --> m1^2/Q^2 = (1-b1) a1  ==> a1(1-b1) =
  //  m1^2/Q^2 p2: (1-a1) p+ +    b1  p- --> m2^2/Q^2 = (1-a1) b1  ==> b1 =
  //  m2^2/Q^2 1/(1-a1)
  //  ==> 0  = a1 [(1-a1)-m2^2/Q^2] - m1^2/Q^2 (1-a1)
  //  ==> 0  = a1 [a1-1+m2^2/Q^2] - m1^2/Q^2 (a1-1) = 0
  //  ==> 0  = Q^2 a1^2 - a1(Q^2+m1^2-m2^2) + m1^2
  //  ==> a1 = 1/(2 Q^2) { Q^2+m1^2-m2^2 + sqrt[(Q^2-m1^2-m2^2)^2+4m2^2m1^2]}
  Vec4D mom1 = p_part1->Momentum(), mom2 = p_part2->Momentum(),
        mom = mom1 + mom2;
  double Q2 = mom.Abs2(), E = sqrt(Q2) / 2.;
  double m1 = p_part1->Flavour().HadMass(), m12 = sqr(m1);
  double m2 = p_part2->Flavour().HadMass(), m22 = sqr(m2);
  double sm = Flavour(newkfc).HadMass();
  if (Q2 < m12 + m22)
    return false;

  ///////////////////////////////////////////////////////////////////
  // 1, Boost to the CoM of the Emitter + Spectator
  // 2, Rotate onto Z
  // 3, Boost to the CoM frame of the emitter only
  // 4, Calculate the momentum of G (M^2-m^2)/(2*m), and sample for cos(theta) and phi.
  // 5, Boost back
  ///////////////////////////////////////////////////////////////////
  Poincare intoCMS = Poincare(mom);
  intoCMS.Boost(mom1);
  intoCMS.Boost(mom2);
  Poincare ontoZ = Poincare(mom1, E * s_AxisP);
  ontoZ.Rotate(mom1);
  ontoZ.Rotate(mom2);
  Poincare intoRESTFRAME = Poincare(mom1);
  intoRESTFRAME.Boost(mom1);

  double z = (m12 - sqr(sm))/(2.*m1);
  if (z<0.) {
  msg_Out()
  	     <<"E1 = "<<mom1[0]<<" - m1 = "<<m1<<" and\n"
  	     <<"m12 = "<<m12<<" Q2 = "<<Q2<<" - m22 = "<<m22<<"\n";
  }
  double p = (m12 - sqr(sm))/(2.*m1);
  double costh = 2.0*ran->Get() - 1.0;
  double sinth = sqrt(1.0 - costh*costh);
  double phi   = 2.0*M_PI*ran->Get();
  Vec3D pg(p*sinth*cos(phi),p*sinth*sin(phi),p*costh);
  m_mom[0][0] = sqrt(sqr(sm)+ p*p);
  m_mom[1][0] = p;
  for(size_t i = 1; i < 4; i++) {
    m_mom[0][i] = -pg[i];
    m_mom[1][i] =  pg[i];
  }
  m_mom[2] = mom2;
  for (size_t i = 0; i < 2; i++) {
    intoRESTFRAME.BoostBack(m_mom[i]);
    ontoZ.RotateBack(m_mom[i]);
    intoCMS.BoostBack(m_mom[i]);
  }
  ontoZ.RotateBack(m_mom[2]);
  intoCMS.BoostBack(m_mom[2]);
  // msg_Out()
  // 	   <<p_part1->Momentum()<<" + "<<p_part2->Momentum()<<"\n"
  // 	   <<m_mom[0]<<" + " <<m_mom[1]<<" + " <<m_mom[2]<<"\n"
  // 	   <<"Check: "<<mom<<" vs. "<<(m_mom[0]+m_mom[1]+m_mom[2])<<"\n";
  return (z>0.);
}
void OctetMeson_Decayer::UpdateColouredObjectsAndAddHadron() {

  kf_code octetkfc = p_part1->Flavour().Kfcode();
  Proto_Particle *meson = new Proto_Particle(Flavour(newkfc), m_mom[0]);
  p_hadrons->push_back(meson);
  p_part1->SetFlavour(Flavour(kf_gluon));
  p_part1->SetMomentum(m_mom[1]);
  p_part2->SetMomentum(m_mom[2]);
  // msg_Out()<<"\n"<<(*meson)<<(*p_part1)<<(*p_part2)<<"\n";
}
