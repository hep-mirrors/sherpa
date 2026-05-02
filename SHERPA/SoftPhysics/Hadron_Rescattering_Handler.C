#include "SHERPA/SoftPhysics/Hadron_Rescattering_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <cmath>
#include <algorithm>

using namespace SHERPA;
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Hadron_Rescattering_Handler::Hadron_Rescattering_Handler() :
  m_name("Hadron_Rescattering"), m_on(true),
  m_rescattering(Hadron_Rescatterings(true))
{
  msg_Out() << "======================================================\n"
            << METHOD << ": Hadron rescattering is "
            << (m_on ? "ON" : "OFF") << "\n"
            << "======================================================\n";
  if (!m_on) return;
  p_coalescence = new Deuteron_Coalescence(
    2,//0.7,                           // b_coal in fm
    0.197326/(2),                          // p_coal in GeV using uncertainty principle : p = hba/(2r)
    CoalescenceModel::CrossSectionModel   // or HardDisk, choose as you wish. 
  );
  
  m_rescattering.Initialize();
  msg_Out() << "Hadron_Rescatterings module initialized!" << std::endl;
  m_particles.clear();
}

Hadron_Rescattering_Handler::~Hadron_Rescattering_Handler() {}

void Hadron_Rescattering_Handler::HarvestParticles(Blob * blob) {
  msg_Out()<<METHOD<<" ----------------------------------------\n";
  if (m_treatedblobs.find(blob)!=m_treatedblobs.end()) {
    msg_Out()<<"* Blob "<<blob->Id()<<" already harvested.\n";
    return;
  }
  for (size_t i=0;i<blob->NOutP();i++) 
  {
    Particle * part1 = blob->OutParticle(i);
    if (m_particles.size()>0) 
    {
      for (std::set<Particle * >::iterator pit=m_particles.begin(); pit!=m_particles.end();pit++) 
      {
        // msg_Out()<<"Scheduling "<<part1->Flav()<<" and "<<(*pit)->Flav()<<std::endl;
	      Particle * part2 = (*pit);
	      Schedule(part1,part2);
      }
    }
    m_particles.insert(part1);
  }

}

void Hadron_Rescattering_Handler::
Schedule(Particle * part1,Particle * part2) {

    //boost to z begin
if (!part1 || !part2) {
    msg_Out() << "Null particle pointer!" << std::endl;
    return;
  }
  
  Vec4D pA = part1->Momentum(), pB = part2->Momentum();
  Vec4D xA = part1->Position(), xB = part2->Position();
  Vec4D P = pA + pB; 
  
  Poincare boost(P);
  boost.Boost(pA); boost.Boost(pB);
  boost.Boost(xA); boost.Boost(xB);
  Vec4D dx = xA-xB;    double b = sqrt(dx[1]*dx[1] + dx[2]*dx[2]);

    // Offset particles to position when the last particle is created.
  double tCreation = std::max(xA[0], xB[0]);
  double zA = xA[3] + (tCreation - xA[0]) * pA[3] / pA[0];
  double zB = xB[3] + (tCreation - xB[0]) * pB[3] / pB[0];
  double tCollision = tCreation + (zB - zA)/(pA[3] / pA[0]-pB[3] / pB[0]);

  // if (zA >= zB) return ;
    //boost to z end

  Vec3D  v1   = Vec3D(part1->Momentum())/part1->Momentum()[0];
  Vec3D  v2   = Vec3D(part2->Momentum())/part2->Momentum()[0];
  double v12  = v1.Sqr(),                 v22 = v2.Sqr();
  double t1   = part1->Position()[0],     t2  = part2->Position()[0];
  Vec3D  x1   = Vec3D(part1->Position()), x2  = Vec3D(part2->Position());

  Vec3D  d12  = x1-t1*v1-x2+t2*v2,        dv = v1-v2;
  double t = -(d12*dv)/(dv.Sqr());
  double D2 = (d12 + t*dv).Sqr(); //mm^2
  Flavour fl1 = part1->Flav(),            fl2 = part2->Flav();

  bool isPair = ((fl1.Kfcode()==2212 && fl2.Kfcode()==2112) ||
                (fl1.Kfcode()==2112 && fl2.Kfcode()==2212));
  if (isPair && ( t>t1 && t>t2)) 
  {
      m_coal_candidates[t] = new Collision(t, D2, part1, part2);
  }
//disregard the collision if D2 > 10 fm squared. 
  if ( t>t1 && t>t2 &&
       (fl1.IsStable()  || (!fl1.IsStable() && t1+rpa->hBar()/fl1.Width()>t)) &&
       (fl2.IsStable() || (!fl2.IsStable() && t2+rpa->hBar()/fl2.Width()>t))) 
  {
    // msg_Out()<<"D2 value is the closest approach distance squared: "<< D2 << std::endl;
    // msg_Out()<<"b is the closest distance: "<< std::pow(b,2)<<std::endl;

    // msg_Out()<<METHOD<<"["
    //        <<std::setw(10)<<part1->Flav()<<": "<<part1->Position()<<",\n"
    //         <<std::string(39,' ')
    //         <<std::setw(10)<<part2->Flav()<<": "<<part2->Position()<<"]:\n"
    //         <<"   - collision time t = "<<t
    // 	      <<", distance^2 = "<<D2 << std::endl;
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
  m_particles.clear();
  m_collisions.clear();

  for (auto &c : m_coal_candidates) delete c.second;
  m_coal_candidates.clear();
}

void Hadron_Rescattering_Handler::Output() {
  msg_Out() << "======================================================\n"
            << "Hadron_Rescattering_Handler::Output(): rescattering is "
            << (m_on ? "ON" : "OFF") << "\n"
            << "======================================================\n";
}

Blob* Hadron_Rescattering_Handler::Rescatter() {
    while (!m_collisions.empty()) {
        Collision* collision = m_collisions.begin()->second;
        m_collisions.erase(m_collisions.begin());

        Particle* pA = collision->p_A;
        Particle* pB = collision->p_B;

        // skip if either particle was already used
        if (pA->Status() != part_status::active) continue;
        if (pB->Status() != part_status::active) continue;

        Blob* blob = m_rescattering(pA, pB, collision->m_dist2);
        msg_Out() << "Handler: collision processed, "
                  << m_collisions.size() << " remaining\n";
        if (blob) return blob;
        delete collision;
    }
    return nullptr;
}

Blob* Hadron_Rescattering_Handler::RescatterCoalescence() {
    while (!m_coal_candidates.empty()) {
      // msg_Out() << "START call: "
      //     << m_coal_candidates.size() << " candidates"<<std::endl;
        Collision* c = m_coal_candidates.begin()->second;
        m_coal_candidates.erase(m_coal_candidates.begin());

        if (c->p_A->Status() != part_status::active) { delete c; continue; }
        if (c->p_B->Status() != part_status::active) { delete c; continue; }

        double b2_fm2 = c->m_dist2 * 1e24;   // mm^2 → fm^2
        Blob* blob = (*p_coalescence)(c->p_A, c->p_B, b2_fm2);
        msg_Out() << "Handler: RescatterCoalescence processed, for particles  "<< c->p_A->Flav() << " and " << c->p_B->Flav()
                   << ", "
                  << m_coal_candidates.size() << " remaining"<<std::endl;

        delete c;

        if (blob) return blob;
    }
    return nullptr;
}


