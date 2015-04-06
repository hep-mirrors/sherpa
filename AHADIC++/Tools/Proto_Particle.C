#include "AHADIC++/Tools/Proto_Particle.H"
#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"

using namespace AHADIC;
using namespace ATOOLS;

namespace AHADIC {
  long int Cluster::s_cluster_count=0;
  long int Cluster::s_cluster_number=0;

  long int control::s_AHAparticles=0;
  long int control::s_AHAprotoparticles=0;
  long int control::s_AHAblobs=0;

  std::list<Proto_Particle *>      Proto_Particle::s_actives;
  std::list<Proto_Particle_List *> Proto_Particle_List::s_actives;
  std::list<ListOfPPLs *>          ListOfPPLs::s_actives;
  std::list<Cluster *>             Cluster::s_actives;
  std::list<Cluster_List *>        Cluster_List::s_actives;
}


Proto_Particle::Proto_Particle(const Proto_Particle & pp) :
  m_flav(pp.m_flav), m_mom(pp.m_mom), m_info(pp.m_info), 
  m_mass(pp.m_mass), m_kt2max(pp.m_kt2max),
  p_partner(pp.p_partner)
{ 
  control::s_AHAprotoparticles++; 
  s_actives.push_back(this);
}

Proto_Particle::
Proto_Particle(Flavour flav,Vec4D mom,char info) :
  m_flav(flav), m_mom(mom), m_info(info), 
  m_mass(hadpars->GetConstituents()->Mass(flav)), m_kt2max(0.), 
  p_partner(NULL)
{ 
  if (!(flav.IsGluon() || flav.IsDiQuark())) {
    if (flav.IsQuark() && flav.Kfcode()>5) {
      std::cerr<<"Error in Proto_Particle::Proto_Particle():\n"
	       <<"   Tried to form a cluster particle from a "<<flav<<".\n"
	       <<"   Please make sure that heavy coloured objects decay "
	       <<"before they enter hadronization.\n"
	       <<"   Will exit the run.\n";
      exit(0);
    }
  }
  control::s_AHAprotoparticles++; 
  s_actives.push_back(this);
}


Proto_Particle::~Proto_Particle()
{ 
#ifdef memchecker
  std::cout<<"### delete Proto_Particle: ("<<m_flav<<"/"<<this<<").\n";
#endif
  control::s_AHAprotoparticles--; 
  s_actives.remove(this);
}

bool Proto_Particle::CheckConsistency(std::ostream & s,std::string method) {
  if (dabs(m_mass-hadpars->GetConstituents()->Mass(m_flav))>1.e-6 ||
      dabs(m_mass-sqrt(m_mom.Abs2()))>1.e-6 ||
      dabs(sqrt(m_mom.Abs2())-hadpars->GetConstituents()->Mass(m_flav))>1.e-6) {
    s<<"Error in "<<METHOD<<" called by "<<method<<":\n"
     <<"   Masses and momenta not consistent for "
     <<m_flav<<"("<<m_mass<<"),"
     <<" sqrt(mom^2) = "<<sqrt(m_mom.Abs2())
     <<" & constituent mass = "<<hadpars->GetConstituents()->Mass(m_flav)
     <<".\n";
    return false;
  }
  return true;
}

std::ostream & AHADIC::
operator<<(std::ostream & s, const Proto_Particle& proto) {
  s<<"   "<<proto.m_info<<" : "<<proto.m_flav<<" "<<proto.m_mom
   <<" "<<sqrt(Max(0.,proto.m_mom.Abs2()))
   <<", kt_max = "<<sqrt(Max(0.,proto.m_kt2max))<<", "
   <<"pt = "<<proto.m_mom.PPerp()<<", y = "<<proto.m_mom.Y()<<std::endl;
  return s;
}

std::ostream & AHADIC::
operator<<(std::ostream & s, const Proto_Particle_List & pl) {
  s<<"Proto_Particle_List with "<<pl.size()<<" elements:\n";
  for (PPL_Const_Iterator pit=pl.begin(); pit!=pl.end(); ++pit) 
    s<<(**pit)<<std::endl;
  return s;
}

