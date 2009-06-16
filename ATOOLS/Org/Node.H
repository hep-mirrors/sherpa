#ifndef ATOOLS_Org_Node_H
#define ATOOLS_Org_Node_H

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

    void operator<<(Node<Node_Type> *const prev);

    Node<Node_Type> *const operator--() const;

  };// end of class Node

}// end of namespace ATOOLS

#include "ATOOLS/Org/Node.C"

#endif
