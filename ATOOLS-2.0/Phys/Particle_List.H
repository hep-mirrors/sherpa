#ifndef Particle_List_H
#define Particle_List_H


#include <deque>
#include <iostream>
// Apacic physics includes
#include "Particle.H"
#include "List_Algorithms.H"

namespace ATOOLS {

  typedef std::deque<Particle*> Particle_List;
  typedef Particle_List::iterator Particle_Iterator;

  std::ostream & operator<<(std::ostream & s, Particle_List & pl);


}

#endif // Particle_List_H
