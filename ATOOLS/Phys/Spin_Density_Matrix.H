#ifndef ATOOLS_Phys_Spin_Density_Matrix_H
#define ATOOLS_Phys_Spin_Density_Matrix_H
#include <iostream>
#include <vector>
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {
  class Spin_Density_Matrix {
  protected:
    //! The SDM is a "m_size x m_size" matrix
    size_t m_size;
    //! Entries of the the matrix
    std::vector<Complex> m_entries;       // entries of the matrix

  public:
    //! Default constructer that does basically nothing
    Spin_Density_Matrix():m_size(0) {};
     
    //! Creates a squared SDM of the form "size x size"
    Spin_Density_Matrix( size_t size );
    //! Creates a squared SDM with the entries <var>entr</var>
    Spin_Density_Matrix( std::vector<Complex> entr );
    //! Destructor
    ~Spin_Density_Matrix() {};

    //! Returns 2 x spin of the matrix: "size-1"
    int     Spin()      { return m_size-1; }        // return 2*spin
     
    //! Returns size of the matrix: "size"
    size_t  Size()      { return m_size; }
     
    //! Returns the number of entries: "size^2"
    size_t  NrEntries() { return m_entries.size(); }
     
    //! Sets the SDM to a trace-one normalised unit matrix
    void    SetNoCorrelation();
     
    //! Sets the SDM to a unit matrix
    void    SetUnitMatrix();
     
    //! Returns the trace of the SDM
    Complex Trace();
     
    //! Normalises the SDM such that it has trace=1
    void    Normalise(); 
     
    //! Used for adding two SDM's
    Spin_Density_Matrix& operator+=(Spin_Density_Matrix SDM);
     
    //! Used for multiplying an SDM with a complex number
    Spin_Density_Matrix& operator*=(Complex factor );
     
    //! Returns a single entry of the SDM (\f$0 \leq entry \leq size^2-1\f$)
    Complex& operator[](size_t entry);
    //! Returns the entry (a,b) of the SDM (\f$0 \leq a,b \leq size-1\f$)
    Complex& operator()(size_t a, size_t b);
    //! Prints the SDM
    void Print(); 
    //! Compares two SDM's
    bool operator== (const Spin_Density_Matrix & sdm );
     
  };

  /*!
    \file Spin_Density_Matrix.H
    \brief Declares the class ATOOLS::Spin_Density_Matrix
    */
  /*!
    \class Spin_Density_Matrix
    \brief Used for spin density matrices ans decay matrices

    The spin density matrices and the decay matrices are  essential tools for the implementation
    of the spin correlation algorithm, see hep-ph/0110108. Both kind of matrices
    are stored as Spin_Density_Matrix. However, any other (squared) matrix can be stored as Spin_Density_Matrix
    as well.
    */
  
}

#endif
