#ifndef ATOOLS_Math_Permutation_H
#define ATOOLS_Math_Permutation_H

namespace ATOOLS {

  class Permutation {
  private:
    int   m_n;
    int** p_an;
    int   *p_per,*p_snum;
  public:
    Permutation(int);                   
    ~Permutation();                   
    int MaxNumber();
    int* Get(int);
  };


}

#endif
