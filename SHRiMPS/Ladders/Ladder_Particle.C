#include "SHRiMPS/Ladders/Ladder_Particle.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Ladder_Particle::Ladder_Particle(const ATOOLS::Flavour & flav,const Vec4D & mom,const Vec4D & pos) :
  m_flav(flav), m_mom(mom), m_pos(pos),
  m_flow(Flow()),
  m_marked(false), m_beam(-1), m_IS(false)
{}

Ladder_Particle::Ladder_Particle(const Ladder_Particle & part) :
  m_flav(part.m_flav), m_mom(part.m_mom), m_pos(part.m_pos),
  m_flow(Flow()),
  m_marked(part.m_marked),m_beam(part.m_beam),m_IS(part.m_IS)
{}

void Ladder_Particle::SetFlow(const unsigned int & pos,const int & code) {
  if(code==-1) m_flow.SetCode(pos);
  else m_flow.SetCode(pos,code);
}


Particle * Ladder_Particle::GetParticle() {
  Particle * part = new Particle(-1,m_flav,m_mom,m_IS?'I':'F');
  part->SetNumber();
  part->SetFlow(1,GetFlow(1));
  part->SetFlow(2,GetFlow(2));
  if (m_beam>=0) part->SetBeam(m_beam);
  return part;
}

namespace SHRIMPS {
std::ostream &
operator<<(std::ostream & s, const Ladder_Particle & part) {
  s<<"   "<<part.Flavour()<<"  "<<part.Momentum()<<" "
   <<"(y="<<part.Momentum().Y()<<", kt^2="<<part.Momentum().PPerp2()<<") "
   <<"{"<<part.GetFlow(1)<<" "<<part.GetFlow(2)<<"}"
   <<" at "<<part.Position()<<".\n";
  return s;
}

std::ostream & operator<<(std::ostream & s, const LadderMap & lmap) {
  size_t i(0);
  s<<"In total "<<lmap.size()<<" emissions:\n";
  for (LadderMap::const_iterator yiter=lmap.begin();yiter!=lmap.end();yiter++) {
    s<<"  y_{"<<(i++)<<"} = "<<yiter->first<<"\n";
  }
  return s;
}
}
