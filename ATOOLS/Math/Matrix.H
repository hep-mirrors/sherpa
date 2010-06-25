#ifndef ATOOLS_Math_Matrix_H
#define ATOOLS_Math_Matrix_H


// (4x4) Matrix class...

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MyComplex.H"

namespace ATOOLS {

  template<int _rank>
  class Matrix {
  protected:
    double ** p_m;
    void NumRecipesNotation();
    void AmegicNotation();
    void Jacobi(double d[], Matrix<_rank>&, int *) const;  
  public:
    Matrix();
    Matrix(const double ma[_rank][_rank]);
    Matrix(const Matrix<_rank>&);
    ~Matrix();
    Matrix<_rank>& operator=(const Matrix<_rank>&);
    Matrix<_rank>  operator*(const double);
    Matrix<_rank>  operator*(const Matrix<_rank>&);
    double       * operator[](int i);
    const double * operator[](int i) const;
    void MatrixOut() const;  
    int Rank() const {return _rank;}
    void Diagonalize(double *,Matrix<_rank>&) const;
    void DiagonalizeSort(double *,Matrix<_rank>&) const;
    Matrix<_rank> Dagger();
    //Vec4D operator*(const Vec4D&); 
  };

  class CMatrix {
  protected:
    Complex ** p_m;
    int        m_rank;
  public:
    CMatrix() { m_rank = 0; p_m = NULL; }
    CMatrix(int);
    CMatrix(Complex**,int);
    CMatrix(const CMatrix&);
    ~CMatrix();
    CMatrix         Conjugate();
    int             Rank() const            { return m_rank; }
    Complex       * operator[](int i)       { return p_m[i]; }
    const Complex * operator[](int i) const { return p_m[i]; }
    Vec4C operator* (const Vec4C& cvec);
    Complex ** GetMatrix()      const { return p_m; } 
  };


  //Global functions.

  template<int _rank>
  inline Matrix<_rank> operator*(const double scal, const Matrix<_rank>& in) { 
    Matrix<_rank> out;
    for(short int i=0; i<_rank; i++) {
      for(short int j=0; j<_rank; j++) {
	out[i][j]=scal*in[i][j];
      }
    }
    return out;
  }		

  template<int _rank>
  Matrix<_rank> operator*(const Matrix<_rank>& a,const Matrix<_rank>& b) {
    Matrix<_rank> out;

    for(short int i=0; i<_rank; i++) {
      for(short int j=0; j<_rank; j++) {
	out[i][j] = 0.;
	for(short int k=0; k<_rank; k++) out[i][j] += a[i][k]*b[k][j];
      }
    }
    return out;
  }

  inline Vec4D operator*(const Matrix<4>& a,const Vec4D& b) {
    Vec4D out;

    for(short int i=0; i<4; i++) {
      out[i] = 0.;
      for(short int k=0; k<4; k++) out[i] += a[i][k]*b[k];
    }
    return out;
  }

  inline Vec4D operator*(const Vec4D& a,const Matrix<4>& b) {
    Vec4D out;

    for(short int i=0; i<4; i++) {
      out[i] = 0.;
      for(short int k=0; k<4; k++) out[i] += a[k]*b[k][i];
    }
    return out;
  }

  CMatrix operator*(const Complex scal, const CMatrix& in);
  CMatrix operator*(const CMatrix& a,const CMatrix& b);

  /*!
    \file 
    \brief contains class Matrix
  */

  /*!
    \class Matrix
    \brief is a template for a \f$n\times n\f$ matrix and some operations on it

    This matrix template can be used to represent any (real) \f$n\times n\f$ 
    (i.e. quadratic) matrix. 
    
    The following operations are provided:
     - Multiplication with a scalar
     - Multiplication with another matrix of same rank
     - Transpose of a matrix
     - Diagonalization of a matrix
     - Diagonalization with elements sorted
     .
     For \f$4\times 4\f$ matrices there is the
     - Mulitplication with a 4-vector (Vec4D)
     .
     also implemented.
  */

  /*!
    \fn void Matrix::NumRecipesNotation();
    \brief changes access mode to matrix
  */

  /*!
    \fn void Matrix::AmegicNotation();
    \brief changes access mode to matrix
  */

  /*!
    \fn void Matrix::Jacobi(double d[], Matrix<_rank>&, int *);  
    \brief uses Jacobi method to obtain eigenvalues and eigenvectors for a matrix
  */

  /*!
    \fn Matrix::Matrix();
    \brief Default constructor, initializes a matrix filed with zeros
  */

  /*!
    \fn Matrix::Matrix(const double ma[_rank][_rank]);
    \brief Special constructor, initializes a matrix with a given C array
  */

  /*!
    \fn Matrix::Matrix(const Matrix<_rank>&);
    \brief Copy constructor
  */

  /*!
    \fn Matrix::~Matrix();
    \brief Destructor
  */

  /*!
    \fn Matrix<_rank>& Matrix::operator=(const Matrix<_rank>&);
    \brief Assignment operator
  */

  /*!
    \fn Matrix<_rank> Matrix::operator*(const double);
    \brief implements the multiplication with a scalar
  */

  /*!
    \fn Matrix<_rank> Matrix::operator*(const Matrix<_rank>&);
    \brief implements the multiplication with a matrix of same rank
  */

  /*!
    \fn double* Matrix::operator[](int i) {return p_m[i];}
    \brief provides write access matrix via p_m[i][j]
  */

  /*!
    \fn const double* Matrix::operator[](int i) const {return p_m[i];}
    \brief provides read only acces to a constant matrix via p_m[i][j]
  */

  /*!
    \fn void Matrix::MatrixOut() const;  
    \brief prints a matrix on screen
  */

  /*!
    \fn const int Matrix::Rank() const {return _rank;}
    \brief returns the rank of the matrix
  */

  /*!
    \fn void Matrix::Diagonalize(double*,Matrix<_rank>&);
    \brief diagonizes the matrix
  */

  /*!
    \fn void Matrix::DiagonalizeSort(double*,Matrix<_rank>&);
    \brief diagonizes the matrix (eigenvalues are sorted)
  */

  /*!
    \fn Matrix<_rank> Matrix::Dagger();
    \brief transposes the matrix
  */

  // other methods

  /*!
    \fn inline Vec4D operator*(const Matrix<4>& a,const Vec4D& b)
    \brief multiplication of a \f$4\times 4\f$ matrix with a 4-vector
  */
  
  /*!
    \fn inline Vec4D operator*(const Vec4D& a,const Matrix<4>& b)
    \brief multiplication of a 4-vector with a \f$4\times 4\f$ matrix
  */
}

#endif

  
