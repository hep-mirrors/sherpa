#include "ATOOLS/Phys/Hadron_Collision.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

Hadron_Collision::Hadron_Collision(Particle * part1,Particle * part2) :
  m_Lorentz(IntoOrientedCoM(part1->Momentum(), part2->Momentum())),
  m_InvLorentz(m_Lorentz), 
  m_on(false)
{
  m_particles[0] = part1;
  m_particles[1] = part2;
  for (size_t i=0;i<2;i++) {
    m_momentaCMS[i]    = m_Lorentz*m_particles[i]->Momentum();
    m_velocitiesCMS[i] = m_momentaCMS[i]/m_momentaCMS[i][0];
    m_positionsCMS[i]  = m_Lorentz*m_particles[i]->Position();
  }
  // Time of closest approach in c.m. system: motion only along the z-axis, hence
  // t_{cms} = [(z1-v1*t1)-(z2-v2*t2)]/(v2-v1)
  // This time has to be larger than the creation times t_1 and t_2 of the two
  // particles to make sure we have the right temporal sequence.
  m_tCMS         = ( ( (m_positionsCMS[1][3]-m_velocitiesCMS[1][3]*m_positionsCMS[1][0]) -
		       (m_positionsCMS[0][3]-m_velocitiesCMS[0][3]*m_positionsCMS[0][0]) )/
		     (m_velocitiesCMS[0][3]-m_velocitiesCMS[1][3]) );
  if (m_tCMS<=m_positionsCMS[1][0] || m_tCMS<=m_positionsCMS[0][0]) return;
  // Rescatter time is smaller than average lifetime of particles if they are unstable.
  // We need to refine this, as particles may well decay before ...
  if ( (!part1->Flav().IsStable() &&
	m_positionsCMS[0][0]+rpa->hBar()/part1->Flav().Width()>m_tCMS) ||
       (!part2->Flav().IsStable() &&
	m_positionsCMS[1][0]+rpa->hBar()/part2->Flav().Width()>m_tCMS) ) return;       
  // In the c.m. system, this minimal distance is given by the impact parameter
  // of the two particles.
  m_b2           = (m_positionsCMS[0]-m_positionsCMS[1]).PPerp2();
  // This gives the "collision" point in the cms as
  // (t_{cms}, (x1+x2)/2, (y1+y2)/2, 0) 
  m_positionCMS  = Vec4D(m_tCMS,(m_positionsCMS[0][1]-m_positionsCMS[1][1])/2.,
			 (m_positionsCMS[0][2]-m_positionsCMS[1][2])/2.,0.);
  // Invert the boost and rotation into the c.m. system ...
  m_InvLorentz.Invert();
  // ... apply it to the position 4-vector, and fix the other characteristics
  m_positionLab   = m_InvLorentz * m_positionCMS;
  m_tLab          = m_positionLab[0];
  m_momentumLab   = part1->Momentum()+part2->Momentum();
  m_s             = m_momentumLab.Abs2();
  m_on            = true;
  //msg_Out()<<METHOD<<"["<<this<<"] for\n  "<<*InPart(0)<<"\n  "<<*InPart(1)<<"\n";
}
