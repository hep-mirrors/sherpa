#ifndef Flavour_List_H
#define Flavour_List_H

// it is probably neccessary to translate this library without the option
// -fno-implicit-templates

// STL includes
//#include <vector>
//#include <list>
#include <deque>

// Apacic physics includes
#include "Flavour.H"

namespace APHYTOOLS {

  //typedef vector<Flavour> Flavour_List;
  //typedef list<Flavour> Flavour_List;
  typedef std::deque<Flavour> Flavour_List;
  typedef Flavour_List::iterator Flavour_Iterator;

  int TestListFunction();
}
#endif // Flavour_List_H

