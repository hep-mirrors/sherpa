#include "Flavour_List.H"
#include "Message.H"

#include <iostream>

using namespace ATOOLS;

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
  
  
  Flavour_List flavs;
  int i;
  flavs.push_back(Flavour(kf_e));
  flavs.push_back(Flavour(kf_e).Bar());
  flavs.push_back(Flavour(kf_d));
  flavs.push_back(Flavour(kf_gluon));
  flavs.push_back(Flavour(kf_d).Bar());
  for (int k=0;k<flavs.size();++k) msg_Out()<<flavs[k]<<std::endl;
  msg_Out()<<"list of "<<flavs.size()<<" Flavours:"<<std::endl;
  for (Flavour_Iterator fliter=flavs.begin();fliter!=flavs.end();++fliter)
    msg_Out()<<*fliter<<std::endl;
  return 0;
};
*/
