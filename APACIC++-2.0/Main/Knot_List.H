#ifndef Knot_List_H
#define Knot_List_H

#include <deque>
// Apacic physics includes
#include "Knot.H"

namespace APACIC {
  typedef std::deque<Knot*> Knot_List;
  typedef Knot_List::iterator Knot_Iterator;
}
#endif
