#include "SHRiMPS/Event_Generation/Single_Diffractive_Event_Generator.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;

Single_Diffractive_Event_Generator::
Single_Diffractive_Event_Generator(Sigma_SD * sigma,Beam_Remnant_Handler * beams,
				   const int & test) :
  p_sigma(sigma), p_beams(beams), 
  m_beam1(ATOOLS::rpa->gen.Beam1()), m_beam2(ATOOLS::rpa->gen.Beam2()),
  m_out1(ATOOLS::Flavour(kf_none)), m_out2(ATOOLS::Flavour(kf_none)), 
  m_p1(ATOOLS::rpa->gen.PBeam(0)), m_p2(ATOOLS::rpa->gen.PBeam(1)),
  m_pl12(Vec3D(m_p1).Sqr()), m_pl22(Vec3D(m_p2).Sqr()),
  m_sign1(double(-1+2*int(m_p1[3]>0))), m_needsboost(false),
  m_accu(1.e-2), m_test(test)
{
  if (test==-1) return;
  // Assume symmetric collisions only
  if ((ATOOLS::Vec3D(m_p1)+ATOOLS::Vec3D(m_p2)).Sqr()>1.e-4) {
    m_needsboost = true;
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   Beamvectors "<<m_p1<<" and "<<m_p2
	       <<" not in c.m. System."<<std::endl
	       <<"   Will terminate the run."<<std::endl;
    exit(1);
  }

  // Assume pp/ppbar collisions only
  if(s_kftable.find(kf_N_1440)==s_kftable.end()) // if not initialised
    s_kftable[kf_N_1440]=new Particle_Info(kf_N_1440_plus,1.44,0.35,0,0,1,1,0,"N(1440)","N(1440)");
  if(s_kftable.find(kf_N_1440_plus)==s_kftable.end()) // if not initialised
    s_kftable[kf_N_1440_plus]=new Particle_Info(kf_N_1440_plus,1.44,0.35,3,0,1,1,0,"N(1440)+","N(1440)+");


  if (m_beam1==ATOOLS::Flavour(kf_p_plus))            
    m_out1 = ATOOLS::Flavour(kf_N_1440_plus);
  else if (m_beam1==ATOOLS::Flavour(kf_p_plus).Bar()) 
    m_out1 = ATOOLS::Flavour(kf_N_1440_plus).Bar();
  if (m_beam2==ATOOLS::Flavour(kf_p_plus))            
    m_out2 = ATOOLS::Flavour(kf_N_1440_plus);
  else if (m_beam2==ATOOLS::Flavour(kf_p_plus).Bar()) 
    m_out2 = ATOOLS::Flavour(kf_N_1440_plus).Bar();

  if (m_out1==ATOOLS::Flavour(kf_none) || 
      m_out2==ATOOLS::Flavour(kf_none)) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   No proton-proton collisions, instead: "<<m_beam1<<" on "<<m_beam2<<"."<<std::endl
	       <<"   Cannot deal with this setup, will terminate the run."<<std::endl;
    exit(1);
  }

//   m_pl12 = m_pl22 = m_pl12+sqr(m_beam1.Mass())-sqr(m_out1.Mass());
}

bool Single_Diffractive_Event_Generator::
SingleDiffractiveEvent(ATOOLS::Blob_List * blobs,const double & xsec) {
  p_beams->InitialiseCollision();
  ATOOLS::Blob * blob(blobs->FindFirst(ATOOLS::btp::Soft_Collision));
  if (blob && blob->Status()==ATOOLS::blob_status::needs_minBias) {
    if (blob->NInP()>0)  {
      msg_Error()<<"Error in "<<METHOD<<": blob has particles."<<std::endl
		 <<(*blob)<<std::endl;
      blob->DeleteInParticles();
    }
    if (blob->NOutP()>0) {
      msg_Error()<<"Error in "<<METHOD<<": blob has particles."<<std::endl
		 <<(*blob)<<std::endl;
      blob->DeleteOutParticles();
    }

    bool mode(FixKinematics());
    Particle * part1in(new Particle(-1,m_beam1,m_p1));
    part1in->SetNumber();
    Particle * part2in(new Particle(-1,m_beam2,m_p2));
    part2in->SetNumber();
    Particle * part1out(new Particle(-1,m_beam1,m_p1out));
    part1out->SetNumber();
    Particle * part2out(new Particle(-1,m_beam2,m_p2out));
    part2out->SetNumber();

    if (mode) {
      part1out->SetFlav(m_out1);
      part1out->SetFinalMass();
    }
    else {
      part2out->SetFlav(m_out2);
      part2out->SetFinalMass();
    }

    blob->AddToInParticles(part1in);
    blob->AddToInParticles(part2in);
    blob->AddToOutParticles(part1out);
    blob->AddToOutParticles(part2out);
    blob->UnsetStatus(ATOOLS::blob_status::needs_minBias);
    blob->SetStatus(ATOOLS::blob_status::needs_hadrondecays);
    blob->SetStatus(ATOOLS::blob_status::needs_beams);
    blob->SetType(ATOOLS::btp::QElastic_Collision);

    return true;
  }
  return false;
}


bool Single_Diffractive_Event_Generator::FixKinematics() {
  bool mode(true);
  double etot(m_p1[0]+m_p2[0]);
  double pt(p_sigma->PT(mode)), pt2(pt*pt);
  double phi(2.*M_PI*ran->Get()), ptx(pt*cos(phi)), pty(pt*sin(phi));
  double e1,e2,pl1,pl2,m1,m2;
  if(mode){
    m1 = m_out1.Mass();
    m2 = m_beam2.Mass();
  }
  else {
    m1 = m_beam1.Mass();
    m2 = m_out2.Mass();
  }
  e1 = (etot*etot+m1*m1-m2*m2)/(2.*etot);
  e2 = etot-e1;
  pl1 = m_sign1*sqrt(e1*e1-pt2-m1*m1);
  pl2 = -m_sign1*sqrt(e2*e2-pt2-m2*m2);
  m_p1out = Vec4D(e1,ptx,pty,pl1);
  m_p2out = Vec4D(e2,-ptx,-pty,pl2);
  return mode;
}
