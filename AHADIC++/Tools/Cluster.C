#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Org/Message.H"

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
