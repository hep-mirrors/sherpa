#include "SHERPA/SoftPhysics/Hadron_Rescattering_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>

using namespace SHERPA;
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Hadron_Rescattering_Handler::Hadron_Rescattering_Handler() :
  m_name("Hadron_Rescattering"), m_on(true), 
  m_rescattering(Hadron_Rescatterings(true))
{
  if (!m_on) return;
  m_rescattering.Initialize();
  m_particles.clear();
}

Hadron_Rescattering_Handler::~Hadron_Rescattering_Handler() {}

void Hadron_Rescattering_Handler::HarvestParticles(Blob * blob) {
  msg_Out()<<METHOD<<" ----------------------------------------\n";
  if (m_treatedblobs.find(blob)!=m_treatedblobs.end()) {
    msg_Out()<<"* Blob "<<blob->Id()<<" already harvested.\n";
    return;
  }
  for (size_t i=0;i<blob->NOutP();i++) {
    Particle * part1 = blob->OutParticle(i);
    if (m_particles.size()>0) {
      for (std::set<Particle * >::iterator pit=m_particles.begin();
	   pit!=m_particles.end();pit++) {
	Particle * part2 = (*pit);
	Schedule(part1,part2);
      }
    } 
    m_particles.insert(part1);
  }
}

void Hadron_Rescattering_Handler::
Schedule(Particle * part1,Particle * part2) {
  Vec3D  v1   = Vec3D(part1->Momentum())/part1->Momentum()[0];
  Vec3D  v2   = Vec3D(part2->Momentum())/part2->Momentum()[0];
  double v12  = v1.Sqr(),                 v22 = v2.Sqr(); 
  double t1   = part1->Position()[0],     t2  = part2->Position()[0];
  Vec3D  x1   = Vec3D(part1->Position()), x2  = part2->Position();
  Vec3D  d12  = x1-t1*v1-x2+t2*v2,        dv = v1-v2;
  double t    = (d12*dv)/(v22-v12);
  Flavour fl1 = part1->Flav(),            fl2 = part2->Flav();
   if (t>t1 && t>t2 &&
       (fl1.IsStable() || (!fl1.IsStable() && t1+rpa->hBar()/fl1.Width()>t)) &&
       (fl2.IsStable() || (!fl2.IsStable() && t2+rpa->hBar()/fl2.Width()>t))) {
    double D2 = (d12+t*dv).Sqr(); 
    msg_Out()<<METHOD<<"["
	     <<std::setw(10)<<part1->Flav()<<": "<<part1->Position()<<",\n"
	     <<std::string(39,' ')
	     <<std::setw(10)<<part2->Flav()<<": "<<part2->Position()<<"]:\n"
	     <<"   - collision time t = "<<t<<", distance^2 = "<<D2<<".\n";
    m_collisions[t] = new Collision(t,D2,part1,part2);
  }
}

bool Hadron_Rescattering_Handler::operator()() {
  while (!m_collisions.empty()) {
    Collision * collision = m_collisions.begin()->second;
    msg_Out()<<METHOD<<"(t = "<<m_collisions.begin()->first<<" vs "
	     <<collision->m_time<<", dist2 = "<<collision->m_dist2<<")\n";
    Blob * blob = m_rescattering(collision->p_A,collision->p_B,
				 collision->m_dist2);
    m_collisions.erase(m_collisions.begin());
    if (blob) break;
  }
  return true;
}

void Hadron_Rescattering_Handler::CleanUp(const size_t & mode) {
  m_treatedblobs.clear();
}

void Hadron_Rescattering_Handler::Output() {}


