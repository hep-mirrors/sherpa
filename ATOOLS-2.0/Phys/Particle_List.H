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
  typedef Particle_List::const_iterator Particle_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Particle_List & pl);

  template <class Predicate> 
  void copy_if(Particle_Iterator first, Particle_Iterator last, std::back_insert_iterator<Particle_List> out, const Predicate & pred) {
    for ( ; first != last; ++first ) { if ( pred(*first) ) out = new Particle(**first); }
  }

}

#endif // Particle_List_H
