#include "Sort.H"
#include "Vector.H"

using namespace ATOOLS;

template <class Type>
inline int Swap(Type & a1, Type & a2) {
  Type h=a1; a1=a2; a2=h;
  return 1;
} 

template <class Type>
void Bubble_Sort::Up(Type * a, int n)  // starting with the biggest
{
  int test=1;
  for (int end=(n-1);((end>0)&&(test));--end) {
    test=0;
    for (int i=0;i<end;++i) {
      if (a[i]<a[i+1]) test+=Swap(a[i],a[i+1]);
    }
  }
}

template <class Type>
void Bubble_Sort::Down(Type * a, int n)
{
  int test=1;
  for (int end=(n-1);((end>0)&&(test));--end) {
    test=0;
    for (int i=0;i<end;++i) {
      if (a[i]>a[i+1]) test+=Swap(a[i],a[i+1]);
    }
  }
}


// specialisation and explicit instatiation
template <>
void Bubble_Sort::Up<Vec4D>(Vec4D * a, int n)  // starting with the biggest
{
  int test=1;
  for (int end=(n-1);((end>0)&&(test));--end) {
    test=0;
    for (int i=0;i<end;++i) {
      if (a[i][0]<a[i+1][0]) test+=Swap(a[i],a[i+1]);
    }
  }
}

template <>
void Bubble_Sort::Down<Vec4D>(Vec4D * a, int n)
{
  int test=1;
  for (int end=(n-1);((end>0)&&(test));--end) {
    test=0;
    for (int i=0;i<end;++i) {
      if (a[i][0]>a[i+1][0]) test+=Swap(a[i],a[i+1]);
    }
  }
}


/*
template <class Type>
int PartitionUp(Type * a, int left, int right) 
{
  // select pivot:
}

template <class Type>
void QuickSort::Up(Type * a, int n, int start=0)  // starting with the biggest
{
}
*/

 
//template void Bubble_Sort::Up<Vec4D>(Vec4D *, int);
template void Bubble_Sort::Up<double>(double *, int);
template void Bubble_Sort::Up<int>(int *, int);
//template void Bubble_Sort::Down<Vec4D>(Vec4D *, int) const ;
template void Bubble_Sort::Down<double>(double *, int);
template void Bubble_Sort::Down<int>(int *, int);
