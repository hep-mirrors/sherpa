#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AHADIC;
using namespace ATOOLS;

std::set<Cluster *> Cluster::s_clusters = std::set<Cluster *>();


Cluster::Cluster(std::pair<Proto_Particle *,Proto_Particle *> parts) :
  m_parts(parts),
  m_momentum(m_parts.first->Momentum()+m_parts.second->Momentum())
{
  s_clusters.insert(this);
}

Cluster::Cluster(Proto_Particle * part1,Proto_Particle * part2) :
  m_momentum(part1->Momentum()+part2->Momentum())
{
  bool barit(false);
  if ((part1->Flavour().IsQuark() && part1->Flavour().IsAnti()) ||
      (part1->Flavour().IsDiQuark() && !part1->Flavour().IsAnti()))
    barit = true;
  m_parts.first  = barit?part2:part1;
  m_parts.second = barit?part1:part2;
  s_clusters.insert(this);
}
    
Cluster::~Cluster() {
  if (s_clusters.find(this)==s_clusters.end()) {
    msg_Error()<<"Did not find cluster ["<<this<<"]\n";
    return;
  }
  s_clusters.erase(this);
}

const ATOOLS::Vec4D Cluster::Velocity() const {
  Vec3D  p3    = Vec3D(m_momentum);
  double mass2 = Max(m_momentum.Abs2(),1.);
  double p32   = p3.Sqr();
  return Vec4D(1.,(p3/sqrt(p32+mass2)));
}

const Vec4D Cluster::Position() const {
  Vec4D pos(0.,0.,0.,0.);
  double E1 = m_parts.first->Momentum()[0], E2 = m_parts.second->Momentum()[0], norm = E1+E2;
  //msg_Out()<<METHOD<<": time difference = "
  //	   <<dabs(m_parts.first->XProd()[0]-m_parts.second->XProd()[0])<<" from\n"
  //	   <<"   "<<std::setw(4)<<m_parts.first->Flavour()<<": "<<m_parts.first->XProd()<<"\n"
  //	   <<"   "<<std::setw(4)<<m_parts.second->Flavour()<<": "<<m_parts.second->XProd()<<"\n"
  //	   <<"   momentum = "<<m_momentum<<", mass = "<<sqrt(m_momentum.Abs2())<<".\n";
  if (dabs(m_parts.first->XProd()[0]-m_parts.second->XProd()[0])<1.e-8) {
    pos = (E1*m_parts.first->XProd()+E2*m_parts.second->XProd())/(E1+E2);
  }
  else if (m_parts.first->XProd()[0]>m_parts.second->XProd()[0]) {
    Vec4D v2   = m_parts.second->Velocity();
    Vec4D pos2 = ( m_parts.second->XProd() +
		   (m_parts.first->XProd()[0]-m_parts.second->XProd()[0]) * v2); 
    pos = (E1*m_parts.first->XProd()+E2*pos2)/(E1+E2);
  }
  else {
    Vec4D v1   = m_parts.first->Velocity();
    Vec4D pos1 = ( m_parts.first->XProd() +
		   (m_parts.second->XProd()[0]-m_parts.first->XProd()[0]) * v1); 
    pos = (E1*pos1+E2*m_parts.second->XProd())/(E1+E2);
  }
  return pos;
}

const double Cluster::DecayTime(const double & mass1,
				const double & mass2) const {
  double M2  = dabs(m_momentum.Abs2()), m12 = sqr(mass1), m22 = sqr(mass2);
  double PS  = M2/sqrt(sqr(M2-m12-m22)+4.*m12*m22);
  double ME2 = m_momentum[0]/Max(dabs(m_momentum.Abs2()),1.);
  return 2.*M_PI*rpa->hBarc() * PS * ME2;
}

const ATOOLS::Vec4D Cluster::DecayPosition(const double & mass1,
					   const double & mass2) const {
  return Position()+DecayTime(mass1,mass2)*Velocity();
}


void Cluster::Clear() {
  delete m_parts.first;
  delete m_parts.second;
}

void Cluster::Reset() {
  if (!s_clusters.empty()) {
    msg_Error()<<METHOD<<" has to erase "<<s_clusters.size()<<" clusters.\n";
    while (!s_clusters.empty()) {
      std::set<Cluster *>::iterator sit = s_clusters.begin();
      s_clusters.erase(sit);
      delete *sit;
    }
  }
}

namespace AHADIC {

  std::ostream& operator<<(std::ostream& str, const Cluster &cluster) {
  str<<"Cluster ["<<cluster.m_parts.first->Flavour()<<", "
     <<cluster.m_parts.second->Flavour()<<"] "
     <<"("<<cluster.m_momentum<<", "
     <<"mass = "<<sqrt(cluster.m_momentum.Abs2())<<", "
     <<"y = "<<cluster.m_momentum.Y()<<")\n";
  return str;
}

std::ostream & operator<<(std::ostream & s, const Cluster_List & cl) {
  Vec4D totmom(0.,0.,0.,0.);
  for (Cluster_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) 
    totmom += (*cit)->Momentum();
  s<<"Cluster List with "<<cl.size()<<" elements, mom = "<<totmom<<":\n";
  for (Cluster_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    s<<**cit<<std::endl;
  }
  return s;
}

}
