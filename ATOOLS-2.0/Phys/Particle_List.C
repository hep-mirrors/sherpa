#include "Particle_List.H"

#include "Particle_Qualifier.H"

using namespace ATOOLS;

template void copy_if<>(Particle_List::iterator, Particle_List::iterator, 
			std::back_insert_iterator<Particle_List>,
			const Is_Gluon &);
template void copy_if<>(Particle_List::iterator, Particle_List::iterator, 
			std::back_insert_iterator<Particle_List>,
			const Is_Photon &);
template void copy_if<>(Particle_List::iterator, Particle_List::iterator, 
			std::back_insert_iterator<Particle_List>,
			const Is_Final_State &);
template void copy_if<>(Particle_List::iterator, Particle_List::iterator, 
			std::back_insert_iterator<Particle_List>,
			const Is_Charged &);

std::ostream &ATOOLS::operator<<(std::ostream &s,const Particle_List &pl) 
{
  s<<"Particle List with "<<pl.size()<<" elements"<<std::endl;
  for (Particle_List::const_iterator pit=pl.begin(); pit!=pl.end(); ++pit) {
    if (*pit!=NULL) s<<**pit<<"\n";
    else s<<"NULL pointer\n";
  }
  return s;
}

Particle_List::Particle_List():
  m_destructor(NULL) {}

Particle_List::Particle_List(const bool destruct):
  m_destructor(destruct?this:NULL) {}

void Particle_List::Clear()
{
  while (!empty()) {
    delete back();
    pop_back();
  }
}
