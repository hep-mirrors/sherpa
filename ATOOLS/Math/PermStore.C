#include "ATOOLS/Math/PermStore.H"

using namespace ATOOLS;



// class PermStore

// constructor
PermStore::PermStore(size_t pnumber) {
  if (pnumber>1) {
    PermStore* nextp;
    for (size_t y=0;y<pnumber;y++) {
      nextp = new PermStore(pnumber-1);
      push_back(nextp);
    }
  }
  else {
    amplitude=Complex(0.,0.);
    colorstr=Complex(0.,0.);
  }
} 

//destructor
PermStore::~PermStore() { 
  if (size()>1) for (size_t y=0;y<size();y++) delete (*this)[y];
}
 
 
void PermStore::PutAmp(size_t* perm,Complex amp) {
  if (size()>1) {
    size_t j=0;
    size_t i;
    for (i=0;perm[i]!=size()-1;i++);
    j=i;
    for (;i<size()-1;i++) perm[i]= perm[i+1]; 
    ((*this)[j])->PutAmp(perm,amp);
  }
  else amplitude=amp;
}

void PermStore::PutColor(size_t* perm,Complex color) {
  if (size()>1) {
    size_t j=0;
    size_t i;
    for (i=0;perm[i]!=size()-1;i++);
    j=i;
    for (;i<size()-1;i++) perm[i]= perm[i+1];
    ((*this)[j])->PutColor(perm,color);
  }
  else colorstr=color;
}


Complex PermStore::GetAmp(size_t* perm) {
  if (size()>1) {
    size_t j=0;
    size_t i;
    for (i=0;perm[i]!=size()-1;i++);
    j=i;
    for (;i<size()-1;i++) perm[i]= perm[i+1]; 
    return ((*this)[j])->GetAmp(perm);
  }
  else return amplitude;
}


Complex PermStore::GetColor(size_t* perm) {
  if (size()>1) {
    size_t j=0;
    size_t i;
    for (i=0;perm[i]!=size()-1;i++);
    j=i;
    for (;i<size()-1;i++) perm[i]= perm[i+1];
    return ((*this)[j])->GetColor(perm);
  }
  else return colorstr;
}






// class PermStoreFast

// constructor
PermStoreFast::PermStoreFast(size_t pnumber,size_t m_pnumber) {
  if (pnumber>1) {
    PermStoreFast* nextp;
    for (size_t y=0;y<m_pnumber;y++) {
      nextp = new PermStoreFast(pnumber-1,m_pnumber);
      push_back(nextp);
    }
  }
  else {
    amplitude=Complex(0.,0.);
    colorstr=Complex(0.,0.);
  }
} 

//destructor
PermStoreFast::~PermStoreFast() { 
  if (size()>1) for (size_t y=0;y<size();y++) delete (*this)[y];
}
 
 
void PermStoreFast::PutAmp(size_t* perm,Complex amp,size_t pnumber) {
  if (pnumber>1) ((*this)[perm[0]])->PutAmp(&perm[1],amp,pnumber-1);
  else amplitude=amp;
}

void PermStoreFast::PutColor(size_t* perm,Complex color,size_t pnumber) {
  if (pnumber>1) ((*this)[perm[0]])->PutColor(&perm[1],color,pnumber-1);
  else colorstr=color;
}


Complex PermStoreFast::GetAmp(size_t* perm,size_t pnumber) {
  if (pnumber>1) return ((*this)[perm[0]])->GetAmp(&perm[1],pnumber-1);
  else return amplitude;
}


Complex PermStoreFast::GetColor(size_t* perm,size_t pnumber) {
  if (pnumber>1) return ((*this)[perm[0]])->GetColor(&perm[1],pnumber-1);
  else return colorstr;
}
