#ifndef ATOOLS_Phys_Particle_List_H
#define ATOOLS_Phys_Particle_List_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/List_Algorithms.H"
#include <deque>
#include <iostream>

namespace ATOOLS {

  class Poincare;

  class Particle_List: public std::deque<Particle*> {
  private:

    class Destructor {
    private:

      Particle_List *p_list;

    public:

      // destructor
      Destructor(): p_list(NULL) {}
      Destructor(Particle_List *const list): p_list(list) {}

      // destructor
      inline ~Destructor() { if (p_list!=NULL) p_list->Clear(); }

    };// end of class Destructor

    Destructor m_destructor;

  public:

    // constructor
    Particle_List();
    Particle_List(const bool destruct);

    // member functions
    void Clear();

    void Boost(Poincare *const boost) const;
    void BoostBack(Poincare *const boost) const;

    void Rotate(Poincare *const rot) const;
    void RotateBack(Poincare *const rot) const;
    
    void Flip() const;

  };// end of class Particle_List

  std::ostream &operator<<(std::ostream &str,const Particle_List &list);

  template <class Predicate> 
  void copy_if(Particle_List::iterator first,Particle_List::iterator last,
	       std::back_insert_iterator<Particle_List> out,
	       const Predicate & pred) 
  { for (;first!=last;++first) if (pred(*first)) out = new Particle(**first); }

}// end of namespace ATOOLS

#endif
