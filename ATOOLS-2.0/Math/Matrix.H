
#ifndef Matrix_h
#define Matrix_h

// (4x4) Matrix class...

#include "Vector.H"

namespace AMATOOLS {

  template<int _rank>
  class Matrix {
  protected:
    double** m;
  public:
    Matrix();
    Matrix(const double ma[_rank][_rank]);
    Matrix(const Matrix<_rank>&);
    ~Matrix();
    Matrix<_rank>& operator=(const Matrix<_rank>&);
    //vec4d operator*(const vec4d&); 
    Matrix<_rank> operator*(const double);
    Matrix<_rank> operator*(const Matrix<_rank>&);
    double* operator[](int i) {return m[i];}
    const double* operator[](int i) const {return m[i];}
    void matrix_out() const;  
    const int rank() const {return _rank;}
    void Num_Recipes_Notation();
    void Amegic_Notation();
    void Diagonalize(double*,Matrix<_rank>&);
    void Diagonalize_Sort(double*,Matrix<_rank>&);
    void jacobi(double d[], Matrix<_rank>&, int *);  
    Matrix<_rank> Dagger();
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

  inline vec4d operator*(const Matrix<4>& a,const vec4d& b) {
    vec4d out;

    for(short int i=0; i<4; i++) {
      out[i] = 0.;
      for(short int k=0; k<4; k++) out[i] += a[i][k]*b[k];
    }
    return out;
  }

  inline vec4d operator*(const vec4d& a,const Matrix<4>& b) {
    vec4d out;

    for(short int i=0; i<4; i++) {
      out[i] = 0.;
      for(short int k=0; k<4; k++) out[i] += a[k]*b[k][i];
    }
    return out;
  }

}

#endif

  
