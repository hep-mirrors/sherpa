#ifndef Sort_H
#define Sort_H

namespace AORGTOOLS {

  class Bubble_Sort {
  public:
    template <class Type>
    static void Up(Type *, int n);  // puts larges element first
    template <class Type>
    static void Down(Type *, int n); // puts smallest element first
  };


  /*!
    \file
    \brief contains the class AORGTOOLS::Bubble_Sort
  */

  /*!
    \class Bubble_Sort
    \brief a simple and rubust sorting algorithm for arrays

    This class provides a simple sorting algoritm for arrays.
    All methods are static, i.e. you can use them without
    creating a object of this class, e.g.:

    \verbatim
      double a[5] = {5., 3., 8., 1., 2.};
      Bubble_Sort::Down(a,5);
    \endverbatim
  */

  /*!
    \fn static void Bubble_Sort::Up(Type *, int n)  
    \brief puts larges element first
  */

  /*!
    \fn static void Bubble_Sort::Down(Type *, int n)
    \brief puts smallest element first
  */

  // fast sorting algorithm, ideal for big arrays (n>15)
  /* not yet implemented
  class QuickSort {
  public:
    template <class Type>
    void Up(Type *, int n, int start=0);
    template <class Type>
    void Down(Type *, int n, int start=0);
  };
  */

}

#endif
