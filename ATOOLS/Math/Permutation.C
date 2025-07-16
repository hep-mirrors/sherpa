#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

Permutation::Permutation(int n) : m_n(n)     
{
  p_per = new int[m_n];
  p_st = new int[m_n];
  m_maxnum=1;
  for(int i=2;i<=m_n;i++) m_maxnum*=i;
}

Permutation::~Permutation()
{
  delete[] p_st;
  delete[] p_per;
}

int* Permutation::Get(int n) 
{
  if (n>m_maxnum) THROW(fatal_error,"Invalid index");
  for(int i=0;i<m_n;++i) {
    p_st[i]=0;
    p_per[i]=i;
  }
  if (n==0) return p_per;
  int i=1, c=0;
  while (i<m_n) {
    if (p_st[i]<i) {
      if (i%2==0) std::swap<int>(p_per[0],p_per[i]);
      else std::swap(p_per[p_st[i]],p_per[i]);
      if (n==++c) return p_per;
      ++p_st[i];
      i=1;
    }
    else {
      p_st[i]=0;
      ++i;
    }
  }
  return p_per;
}
