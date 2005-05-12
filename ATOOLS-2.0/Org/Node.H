#ifndef Node_H
#define Node_H

#include <vector>

namespace ATOOLS {

  template <class Node_Type>
  class Node: public std::vector<Node_Type> {
  private:
    
    std::vector<Node<Node_Type>*> *p_next;

    Node<Node_Type> *p_previous; 

    Node();
    
  public:

    // constructor
    Node(const Node_Type &node,const bool create=false);
    
    // destructor
    ~Node();
    
    // member functions
    std::vector<Node<Node_Type>*> *operator->();
    std::vector<Node<Node_Type>*> &operator()();

    // inline functions
    inline void operator<<(Node<Node_Type> *const prev) { p_previous=prev; }

    inline Node<Node_Type> *const operator--() const { return p_previous; }

  };// end of class Node

}// end of namespace ATOOLS

#include "Node.C"

#endif
