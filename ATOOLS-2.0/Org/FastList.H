#ifndef FastLists_H
#define FastLists_H

#include <iostream>
#include "Run_Parameter.H"
#include "Message.H"



namespace ATOOLS {
  
  /*
  template<class>
  class ListIterator;
  */  

  template<class F>
  class FastList {
    struct Module {
      Module() : next(0), datapointer(0) {};
      Module(F * f): next(0), datapointer(f) {};
      long int no;
      Module * next;
      F * datapointer;
    };
    Module * root, * last;
    long counter;
  public:
    FastList() : root(0), last(0), counter(0) {};
    ~FastList() {
      Module * iter = root;
      Module * iter1;
      while (iter) {
	iter1 = iter->next;
	delete iter->datapointer;
	delete iter;
	iter = iter1;
      }
    }
    long GetLength() {return counter;};
    void Append(F * f) {
      if (last) {
	last->next = new Module(f);
	last = last->next;
	last->no=++counter;
      } 
      else {
	root = new Module(f);
	last = root;
	last->no=++counter;
      };
    };
    F * Get(long int i) {
      if ((i<=counter)&&(i>0)) {
	Module * act=root;
	while (act->no!=i) {act=act->next;}
	return act->datapointer;
      } 
      msg.Error()<<"Error in FastList: Module no out of bounds!"<<std::endl;
      return 0;
    };
    F * operator[](long int i) { return Get(i+1);};
    friend class Iterator;
    class Iterator {
      Module * root, * current;
    public:
      Iterator(): root(0), current(0) {};
      Iterator(FastList & list): root(list.root), current(list.root) {};
      void Attach(FastList & list) {
	root=list.root;
	current=root;
      };
      F * operator()() { 
	if (current) return current->datapointer;
	else return 0;
      };
      void operator++() { // refix operator
	if (current) current=current->next;
      };
    };
  };


  template<class F>  std::ostream & operator<< (std::ostream & s, FastList<F> & list);
}
#endif
