#include "Blob_List.H"


using namespace APHYTOOLS;

//explicite template instatiatation
template class std::deque<Blob*>;
template class std::back_insert_iterator<Blob_List>;

std::ostream & APHYTOOLS::operator<<(std::ostream & s,const Blob_List & bl) {
  s<<"Blob List with "<<bl.size()<<" elements"<<std::endl;
  for (Blob_Const_Iterator bit=bl.begin(); bit!=bl.end(); ++bit) {
    s<<(*bit)<<std::endl; // note still pointer of blob!
  }
  return s;
}
 
