#include "Particle_List.H"
//#include <iostream>

using namespace ATOOLS;

//explicite template instatiatation
template class std::deque<Particle*>;
template class std::back_insert_iterator<Particle_List>;

namespace ATOOLS {
  template void copy_if<>( Particle_Iterator, Particle_Iterator , 
				      std::back_insert_iterator<Particle_List> ,Is_Gluon );
  template void copy_if<>( Particle_Iterator, Particle_Iterator , 
				      std::back_insert_iterator<Particle_List> ,Is_Photon);
  template void copy_if<>( Particle_Iterator, Particle_Iterator , 
				      std::back_insert_iterator<Particle_List> ,Is_Final_State);
  template void copy_if<>( Particle_Iterator, Particle_Iterator , 
				      std::back_insert_iterator<Particle_List> ,Is_Charged);

  std::ostream & operator<<(std::ostream & s, const Particle_List & pl) {
    s<<"Particle List with "<<pl.size()<<" elements"<<std::endl;
    for (Particle_List::const_iterator pit=pl.begin(); pit!=pl.end(); ++pit) {
      s<<**pit<<std::endl;
    }
    return s;
  }

}
 
