#ifndef Node_H
#define Node_H

#include <vector>

namespace ATOOLS {

  template <class Node_Type>
  class Node: public std::vector<Node_Type> {
  private:
    
    std::vector<Node<Node_Type>*> *p_next;
    
  public:

    // constructor
    Node(const Node_Type &node,const bool create=false);
    
    // destructor
    ~Node();
    
    // member function
    std::vector<Node<Node_Type>*> *operator->();

  };// end of class Node

}// end of namespace ATOOLS

#include "Node.C"

#endif
