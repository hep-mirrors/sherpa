#ifndef Sort_H
#define Sort_H

namespace AORGTOOLS {

  //! a simple and rubust sorting algorithm
  /*! 
    all methods are static, i.e. you can use them without
    creating a object of this class 
  */
  class BubbleSort {
  public:
    template <class Type>
    static void Up(Type *, int n);  // puts larges element first
    template <class Type>
    static void Down(Type *, int n); // puts smallest element first
  };


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
