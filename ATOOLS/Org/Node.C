#ifndef Node_C
#define Node_C

namespace ATOOLS {

  template <class Node_Type>
  Node<Node_Type>::Node(const Node_Type &node,const bool create): 
    std::vector<Node_Type>(1,node), 
    p_next(NULL),
    p_previous(NULL)
  {
    if (create) p_next = new std::vector<Node<Node_Type>*>();
  }
  
  template <class Node_Type>
  Node<Node_Type>::~Node() 
  {
    if (p_next!=NULL) {
      for (typename std::vector<Node<Node_Type>*>::iterator nit=p_next->begin();
	   nit!=p_next->end();++nit) {
	delete (*nit);
      }
      delete p_next;
    }
  }
  
  template <class Node_Type>
  std::vector<Node<Node_Type>*> *Node<Node_Type>::operator->() 
  { 
    return p_next; 
  }
  
  template <class Node_Type>
  std::vector<Node<Node_Type>*> &Node<Node_Type>::operator()() 
  { 
    return *p_next; 
  }
  
  template <class Node_Type>
  void Node<Node_Type>::operator<<(Node<Node_Type> *const prev) 
  { 
    p_previous=prev; 
  }

  template <class Node_Type>
  Node<Node_Type>* Node<Node_Type>::operator--() const 
  { 
    return p_previous; 
  }

}// end of namespace ATOOLS

#endif
