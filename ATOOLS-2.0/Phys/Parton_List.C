#include "Parton_List.H"
//#include <iostream>

using namespace ATOOLS;

//explicite template instatiatation
template class std::deque<Parton*>;
template class std::back_insert_iterator<Parton_List>;

namespace ATOOLS {
  template void copy_if<>( Parton_Iterator, Parton_Iterator , 
				      std::back_insert_iterator<Parton_List> ,Is_Gluon );
  template void copy_if<>( Parton_Iterator, Parton_Iterator , 
				      std::back_insert_iterator<Parton_List> ,Is_Photon);
  template void copy_if<>( Parton_Iterator, Parton_Iterator , 
				      std::back_insert_iterator<Parton_List> ,Is_Final_State);
}

std::ostream & operator<<(std::ostream & s, Parton_List & pl) {
  s<<"Parton List with "<<pl.size()<<" elements"<<std::endl;
  for (Parton_Iterator pit=pl.begin(); pit!=pl.end(); ++pit) {
    s<<(*pit)<<std::endl; // note still pointer of parton!
  }
  return s;
}
 
