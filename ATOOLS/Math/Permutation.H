#ifndef ATOOLS_Math_Permutation_H
#define ATOOLS_Math_Permutation_H

namespace ATOOLS {

  class Permutation {
  private:
    int m_n, *p_per, *p_st;
    long int m_maxnum;
  public:
    Permutation(int);                   
    ~Permutation();                   
    int* Get(int);
    inline int MaxNumber() const { return m_maxnum; }
  };


}

#endif
