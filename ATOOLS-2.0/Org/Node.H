#ifndef Node_H
#define Node_H

#include <vector>

namespace ATOOLS {

  template <class Node_Type>
  class Node: public std::vector<Node_Type> {
  private:
    
    std::vector<Node<Node_Type>*> *p_next;
    
  public:

    Node(const Node_Type &node,const bool create=false): 
      std::vector<Node_Type>(1,node), 
      p_next(NULL) 
    {
      if (create) p_next = new std::vector<Node<Node_Type>*>();
    }
    
    ~Node() 
    {
      if (p_next!=NULL) {
	for (typename std::vector<Node<Node_Type>*>::iterator nit=p_next->begin();
	     nit!=p_next->end();++nit) {
	  delete (*nit);
	}
	delete p_next;
      }
    }
    
    inline std::vector<Node<Node_Type>*> *operator->() { return p_next; }

  };// end of class Node

}// end of namespace ATOOLS

#endif
