#include "Flavour_List.H"
#include "Message.H"

#include <iostream>


using namespace APHYTOOLS;
using namespace AORGTOOLS;

//explicite template instatiatations
//template class list<Flavour>;
template class std::deque<Flavour>;

/*

 // a demo how to use the Flavour iterator :

#include <vector>

//test function for list operations:
int TestListFunction() {
  // make list tests:
  
  std::vector<int> v1;
  v1.push_back(2);
  v1.push_back(3);
  v1.push_back(1);
  
  msg.Out()<<" std::vector<int> :";
  for (int k=0;k<v1.size();++k) msg.Out()<<v1[k]<<" ";
  msg.Out()<<std::endl;
  
  
  Flavour_List flavs;
  int i;
  flavs.push_back(Flavour(kf::e));
  flavs.push_back(Flavour(kf::e).bar());
  flavs.push_back(Flavour(kf::d));
  flavs.push_back(Flavour(kf::gluon));
  flavs.push_back(Flavour(kf::d).bar());
  for (int k=0;k<flavs.size();++k) msg.Out()<<flavs[k]<<std::endl;
  msg.Out()<<"list of "<<flavs.size()<<" Flavours:"<<std::endl;
  for (Flavour_Iterator fliter=flavs.begin();fliter!=flavs.end();++fliter)
    msg.Out()<<*fliter<<std::endl;
  return 0;
};
*/
