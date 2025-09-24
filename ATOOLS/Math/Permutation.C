#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

Permutation::Permutation(int n,int algo):
  m_algo(algo), m_n(n)
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

void Permutation::Swap(int i,int j)
{
  std::swap<int>(p_per[i],p_per[j]);
  std::swap<int>(p_st[i],p_st[j]);
  if ((i==0 && p_st[i]<0) || (i==m_n-1 && p_st[i]>0)) p_st[i]=0;
  if ((j==0 && p_st[j]<0) || (j==m_n-1 && p_st[j]>0)) p_st[j]=0;
}

int Permutation::LargestMobile()
{
  int l(0), lv(-1);
  for (size_t i(0);i<m_n;++i)
    if (p_st[i] && p_per[i]>lv) lv=p_per[l=i];
  return l;
}

int* Permutation::Get(int n) 
{
  if (n>m_maxnum) THROW(fatal_error,"Invalid index");
  if (m_algo==1) {
    for(size_t i=0;i<m_n;++i) {
      p_per[i]=i;
      p_st[i]=-1;
    }
    p_st[0]=0;
    int count(0), m, l;
    if (count==n) return p_per;
    while (p_st[m=LargestMobile()]!=0) {
      Swap(m,l=m+p_st[m]);
      for(size_t i(0);i<m_n;++i)
	if(p_per[i]>p_per[l]) p_st[i]=i<l?+1:-1;
      if (++count==n) return p_per;
    }
    return NULL;
  }
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
