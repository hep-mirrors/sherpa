#include "SINIC++/Event_Generation/Elastic_Event_Generator.H"
#include "SINIC++/Tools/Special_Functions.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace SINIC;

Elastic_Event_Generator::
Elastic_Event_Generator(Sigma_Elastic * sigma,Beam_Remnant_Handler * beams,
			const int & test) :
  p_sigma(sigma), p_beams(beams),
  m_beam1(ATOOLS::rpa->gen.Beam1()), m_beam2(ATOOLS::rpa->gen.Beam2()),
  m_p1(ATOOLS::rpa->gen.PBeam(0)), m_p2(ATOOLS::rpa->gen.PBeam(1)),
  m_pl12(Vec3D(m_p1).Sqr()), m_pl22(Vec3D(m_p2).Sqr()),
  m_sign1(double(-1+2*int(m_p1[3]>0))), m_needsboost(false),
  m_accu(1.e-2), m_test(test)
{
  if (test==-1) return;
  if ((Vec3D(m_p1)+Vec3D(m_p2)).Sqr()>1.e-4) {
    m_needsboost = true;
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   Beamvectors "<<m_p1<<" and "<<m_p2
	       <<" not in c.m. System."<<std::endl
	       <<"   Will terminate the run."<<std::endl;
    exit(1);
  }
}

bool Elastic_Event_Generator::
ElasticEvent(ATOOLS::Blob_List * blobs,const double & xsec) {
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

    FixKinematics();
    Particle * part1in(new Particle(-1,m_beam1,m_p1));
    part1in->SetNumber();
    Particle * part2in(new Particle(-1,m_beam2,m_p2));
    part2in->SetNumber();
    Particle * part1out(new Particle(-1,m_beam1,m_p1out));
    part1out->SetNumber();
    Particle * part2out(new Particle(-1,m_beam2,m_p2out));
    part2out->SetNumber();

    blob->AddToInParticles(part1in);
    blob->AddToInParticles(part2in);
    blob->AddToOutParticles(part1out);
    blob->AddToOutParticles(part2out);
    blob->UnsetStatus(ATOOLS::blob_status::needs_minBias);
    blob->SetStatus(ATOOLS::blob_status::needs_beams);
    blob->SetType(ATOOLS::btp::QElastic_Collision);

    return true;
  }
  return false;
}


void Elastic_Event_Generator::FixKinematics() {
  double pt(p_sigma->PT()), pt2(pt*pt);
  double phi(2.*M_PI*ran->Get()), ptx(pt*cos(phi)), pty(pt*sin(phi));
  double pl1(m_sign1*sqrt(m_pl12-pt2)), pl2(-m_sign1*sqrt(m_pl22-pt2));

  m_p1out = Vec4D(m_p1[0],ptx,pty,pl1);
  m_p2out = Vec4D(m_p2[0],-ptx,-pty,pl2);
}
