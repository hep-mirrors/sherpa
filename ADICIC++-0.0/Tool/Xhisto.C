//bof
//Version: 3 ADICIC++-0.0/2005/07/28

//Implementation of Xhisto.H.



#include <stdio.h>
#include <cassert>
#include "Xhisto.H"





using namespace std;
using namespace ADICIC;





//=============================================================================



Xhisto::Xhisto(size_t stop)
  : m_i(0), m_stop(stop), l_sets() {}



void Xhisto::Reset() {
  m_i=0;
  l_sets.clear();
}



void Xhisto::Insert(const Multidouble& m) {
  //l_sets.size() is a slow operation!
  if(m_i<m_stop) { l_sets.push_back(m); ++m_i;}
}



void Xhisto::Output(const std::string& name) const {
  std::ofstream ofile;
  ofile.open(name.c_str());
  for(std::list<Multidouble>::const_iterator it=l_sets.begin();
      it!=l_sets.end(); ++it) {
    for(size_t i=0; i<it->size(); ++i) ofile<<(*it)[i]<<"  ";
    ofile<<std::endl;
  }
  ofile.close();
}



const size_t Xhisto::Entries() const {
  assert(m_i==l_sets.size());
  return l_sets.size();
}



//=============================================================================





//eof
