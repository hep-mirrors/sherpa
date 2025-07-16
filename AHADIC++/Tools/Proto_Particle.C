#include "AHADIC++/Tools/Proto_Particle.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

std::set<Proto_Particle *> Proto_Particle::s_protos =
  std::set<Proto_Particle *>();

Proto_Particle::Proto_Particle(const Proto_Particle & proto) :
  m_flav(proto.m_flav), m_momentum(proto.m_momentum), m_gen(1),
  m_kt2max(proto.KT2_Max()),
  m_isleading(proto.m_isleading), m_isbeam(proto.m_isbeam)
{
  s_protos.insert(this);
}

Proto_Particle::Proto_Particle(const ATOOLS::Particle & part) :
  m_flav(part.Flav()), m_momentum(part.Momentum()), m_kt2max(sqr(m_momentum[0])),
  m_isleading(false), m_isbeam(part.Info()=='B')
{
  s_protos.insert(this);
}

Proto_Particle::Proto_Particle(const ATOOLS::Flavour & flav,
			       const ATOOLS::Vec4D & mom,
			       bool leading, bool beam) :
  m_flav(flav), m_momentum(mom),m_isleading(leading),m_isbeam(beam)
{
  s_protos.insert(this);
}
  
Proto_Particle::~Proto_Particle()
{
  if (s_protos.find(this)==s_protos.end()) {
    msg_Error()<<"Did not find Proto_Particle ["<<this<<"]\n";
    return;
  }
  s_protos.erase(this);
} 

void Proto_Particle::Reset() {
  for (auto it = s_protos.begin(); it != s_protos.end(); ) {
    delete (*(it++));
  }
}

Particle * Proto_Particle::operator()() {
  return new Particle(-1,m_flav,m_momentum,'P');
}


std::ostream& AHADIC::operator<<(std::ostream & str,
				 const Proto_Particle & proto) {
  str<<"Proto_Particle ["<<proto.Flavour()<<"] "
     <<"("<<proto.Momentum()<<", "
     <<"mass = "<<sqrt(proto.Momentum().Abs2())<<", "
     <<"y = "<<proto.Momentum().Y()<<")\n";
  return str;
}
