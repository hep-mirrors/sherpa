#include "FastList.H"

template<class F>
std::ostream & ATOOLS::operator<< (std::ostream & s, ATOOLS::FastList<F> & list) {
  for (ATOOLS::FastList<F>::Iterator iter(list); iter(); ++iter){
    s<<(*iter())<<endl;
  }  
  return s;
};



// explicit instantiation

//template class FastList<double>;
//template std::ostream & operator<<<double> (std::ostream &, FastList<double> &);
