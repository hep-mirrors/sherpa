#include "AHADIC++/Tools/Proto_Particle.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Proto_Particle::Proto_Particle(const Proto_Particle & proto) :
  m_flav(proto.m_flav), m_momentum(proto.m_momentum),
  m_isleading(proto.m_isleading), m_isbeam(proto.m_isbeam) {
  //msg_Out()<<METHOD<<"(1:"<<this<<", "<<m_flav<<").\n";
}

Proto_Particle::Proto_Particle(const ATOOLS::Particle & part) :
  m_flav(part.Flav()), m_momentum(part.Momentum()),
  m_isleading(false), m_isbeam(part.Info()=='B') {
  //msg_Out()<<METHOD<<"(2:"<<this<<", "<<m_flav<<").\n";
}


Proto_Particle::Proto_Particle(const ATOOLS::Flavour & flav,
			       const ATOOLS::Vec4D & mom,
			       bool leading, bool beam) :
  m_flav(flav), m_momentum(mom),m_isleading(leading),m_isbeam(beam) {
  //msg_Out()<<METHOD<<"(3:"<<this<<", "<<m_flav<<").\n";
}
  
Proto_Particle::~Proto_Particle() {
  //msg_Out()<<METHOD<<"("<<this<<", "<<m_flav<<").\n";
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
