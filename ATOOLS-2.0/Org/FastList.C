#include "FastList.H"

template<class F>
std::ostream & AORGTOOLS::operator<< (std::ostream & s, AORGTOOLS::FastList<F> & list) {
  for (AORGTOOLS::FastList<F>::Iterator iter(list); iter(); ++iter){
    s<<(*iter())<<endl;
  }  
  return s;
};



// explicit instantiation

//template class AORGTOOLS::FastList<double>;
//template std::ostream & AORGTOOLS::operator<<<double> (std::ostream &, AORGTOOLS::FastList<double> &);
