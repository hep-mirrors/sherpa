#include "AHADIC++/Tools/Proto_Particle.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

std::set<Proto_Particle *> Proto_Particle::s_protos =
  std::set<Proto_Particle *>();

Proto_Particle::Proto_Particle(const Proto_Particle & proto) :
  m_flav(proto.m_flav), m_momentum(proto.m_momentum),
  m_xprod(proto.m_xprod), m_xdec(proto.m_xdec), 
  m_gen(1),
  m_kt2max(proto.KT2_Max()),
  m_isleading(proto.m_isleading), m_isbeam(proto.m_isbeam)
{
  s_protos.insert(this);
}

Proto_Particle::Proto_Particle(const ATOOLS::Particle & part) :
  m_flav(part.Flav()), m_momentum(part.Momentum()),
  m_xprod(part.Position()), m_xdec(Vec4D()), 
  m_kt2max(sqr(m_momentum[0])),
  m_isleading(false), m_isbeam(part.Info()=='B')
{
  s_protos.insert(this);
}

Proto_Particle::Proto_Particle(const ATOOLS::Flavour & flav,
			       const ATOOLS::Vec4D & mom,
			       bool leading, bool beam) :
  m_flav(flav), m_momentum(mom),
  m_xprod(Vec4D()), m_xdec(Vec4D()), 
  m_isleading(leading),m_isbeam(beam)
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
  Particle * part = new Particle(-1,m_flav,m_momentum,'P');
  part->SetPosition(m_xprod);
  msg_Out()<<METHOD<<" ["<<std::setw(12)<<part->Flav()<<", "
	   <<std::setw(12)<<(part->Momentum()[0]/part->Flav().HadMass())<<"]: "
	   <<part->Position()<<"\n";
  return part;
}

const ATOOLS::Vec4D Proto_Particle::Velocity() const {
  Vec3D p3     = Vec3D(m_momentum);
  double mass2 = Max(sqr(hadpars->GetConstituents()->Mass(m_flav)),1.);
  double p32   = p3.Sqr();
  return Vec4D(1.,(p3/sqrt(p32+mass2)));
}


std::ostream& AHADIC::operator<<(std::ostream & str,
				 const Proto_Particle & proto) {
  str<<"Proto_Particle ["<<proto.Flavour()<<"] "
     <<"("<<proto.Momentum()<<", "
     <<"mass = "<<sqrt(proto.Momentum().Abs2())<<", "
     <<"y = "<<proto.Momentum().Y()<<")\n";
  return str;
}
