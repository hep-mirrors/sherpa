#ifndef Spin_Density_Matrix_H
#define Spin_Density_Matrix_H
#include <iostream>
#include <vector>
#include <MyComplex.H>
#include <Vector.H>
#include "Flavour.H"

namespace ATOOLS {
  class Spin_Density_Matrix {
  protected:
    //! SDM is a m_size * m_size matrix
    size_t m_size;
    std::vector<Complex> m_entries;       // entries of the matrix

  public:
    Spin_Density_Matrix():m_size(0) {};
    Spin_Density_Matrix( size_t size );
    Spin_Density_Matrix( std::vector<Complex> entr );
    //! Creates a square-sdm with size_a*size_b = total_size
    //Spin_Density_Matrix( size_t total_size );
    ~Spin_Density_Matrix() {};

    int     Spin()      { return m_size-1; }        // return 2*spin
    size_t  Size()      { return m_size; }
    size_t  NrEntries() { return m_entries.size(); }
    void    SetNoCorrelation();
    void    SetUnitMatrix();
    Complex Trace();
    void    Normalise(); 
    double  Contract( std::vector<Complex> * ampls, std::vector<int> * ind ); 
     
    Spin_Density_Matrix& operator+=(Spin_Density_Matrix SDM);
    Spin_Density_Matrix& operator*=(Complex factor );
    Complex& operator[](size_t entry);
    Complex& operator()(size_t a, size_t b);
    void Print(); 
    bool operator== (const Spin_Density_Matrix & sdm );
    friend std::ostream &operator<<
      (std::ostream &ostr, Spin_Density_Matrix &sdm);
     
  };
  
  std::ostream &operator<<(std::ostream &ostr,
      Spin_Density_Matrix &sdm);
  
}

#endif
