#ifndef Blob_List_H
#define Blob_List_H

#include <deque>
#include <iostream>
#include "Blob.H"
#include "List_Algorithms.H"

namespace ATOOLS {
  typedef std::deque<Blob*> Blob_List;
  typedef Blob_List::iterator Blob_Iterator;
  typedef Blob_List::const_iterator Blob_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Blob_List & bl);
}

#endif // Blob_List_H
