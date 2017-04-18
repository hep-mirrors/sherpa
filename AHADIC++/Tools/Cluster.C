#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

Cluster::Cluster(std::pair<Proto_Particle *,Proto_Particle *> parts) :
  m_parts(parts),
  m_momentum(m_parts.first->Momentum()+m_parts.second->Momentum())
{
  //msg_Out()<<METHOD<<"(1:"<<this<<") "
  //	   <<"["<<m_parts.first->Flavour()<<", "
  //	   <<m_parts.second->Flavour()<<"].\n";
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
  //msg_Out()<<METHOD<<"(2:"<<this<<") "
  //	   <<"["<<m_parts.first->Flavour()<<", "
  //	   <<m_parts.second->Flavour()<<"].\n";
}
    
Cluster::~Cluster() {
  //msg_Out()<<METHOD<<"(:"<<this<<").\n";
}

void Cluster::Clear() {
  //msg_Out()<<METHOD<<"(:"<<this<<") "
  //	   <<"["<<m_parts.first->Flavour()<<", "
  //	   <<m_parts.second->Flavour()<<"].\n";
  delete m_parts.first;
  delete m_parts.second;
}


std::ostream& AHADIC::operator<<(std::ostream& str, const Cluster &cluster) {
  str<<"Cluster ["<<cluster.m_parts.first->Flavour()<<", "
     <<cluster.m_parts.second->Flavour()<<"] "
     <<"("<<cluster.m_momentum<<", "
     <<"mass = "<<sqrt(cluster.m_momentum.Abs2())<<", "
     <<"y = "<<cluster.m_momentum.Y()<<")\n";
  return str;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const Cluster_List & cl) {
  Vec4D totmom(0.,0.,0.,0.);
  for (Cluster_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) 
    totmom += (*cit)->Momentum();
  s<<"Cluster List with "<<cl.size()<<" elements, mom = "<<totmom<<":\n";
  for (Cluster_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    s<<**cit<<std::endl;
  }
  return s;
}

