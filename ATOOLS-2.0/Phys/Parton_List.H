#ifndef Parton_List_H
#define Parton_List_H


#include <deque>
#include <iostream>
// Apacic physics includes
#include "Parton.H"
#include "List_Algorithms.H"

namespace ATOOLS {

  typedef std::deque<Parton*> Parton_List;
  typedef Parton_List::iterator Parton_Iterator;

  std::ostream & operator<<(std::ostream & s, Parton_List & pl);


}

#endif // Parton_List_H
