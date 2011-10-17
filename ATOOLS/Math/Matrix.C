#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>


namespace ATOOLS {
  CMatrix operator*(const Complex scal, const CMatrix& in) { 
    CMatrix out(in.Rank());
    for(short int i=0; i<in.Rank(); i++) {
      for(short int j=0; j<in.Rank(); j++) {
	out[i][j]=scal*in[i][j];
      }
    }
    return out;
  }		
  
  CMatrix operator*(const CMatrix& a,const CMatrix& b) {
    if (a.Rank()!=b.Rank()) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"  Tried to multiply two matrices of different rank."<<std::endl
		 <<"  Return 0 and hope for the best."<<std::endl;
      return CMatrix(0);
    }
    CMatrix out(a.Rank());
    for(short int i=0; i<a.Rank(); i++) {
      for(short int j=0; j<a.Rank(); j++) {
	out[i][j] = 0.;
	for(short int k=0; k<a.Rank(); k++) out[i][j] += a[i][k]*b[k][j];
      }
    }
    return out;
  }
}

using namespace ATOOLS;
using namespace std;

template<int _rank>
Matrix<_rank>::Matrix() 
{
  p_m = new double*[_rank];
  for(short int i=0; i<_rank; i++) {
    p_m[i] = new double[_rank];
    for(short int j=0; j<_rank; j++) 
      p_m[i][j]=0.0;
  }
}

template<int _rank>
Matrix<_rank>::Matrix(const double ma[_rank][_rank]) 
{
  p_m = new double*[_rank];
  for(short int i=0; i<_rank; i++) {
    p_m[i] = new double[_rank];
    for(short int j=0; j<_rank; j++) 
      p_m[i][j]=ma[i][j];
  }
}

template<int _rank>
Matrix<_rank>::Matrix(const Matrix<_rank>& in) 
{
  p_m = new double*[_rank];
  for(short int i=0; i<_rank; i++) {
    p_m[i] = new double[_rank];
    for(short int j=0; j<_rank; j++) 
      p_m[i][j]=in[i][j];
  }
}

template<int _rank>
Matrix<_rank>::~Matrix() 
{
  for(short int i=0; i<_rank; i++) delete[] p_m[i];
  delete[] p_m;
}

template<int _rank>
double *Matrix<_rank>::operator[](int i)
{
  return p_m[i]; 
}

template<int _rank>
const double *Matrix<_rank>::operator[](int i) const 
{
  return p_m[i]; 
}

template<int _rank>
Matrix<_rank>& Matrix<_rank>::operator=(const Matrix<_rank>& in) 
{
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) 
      p_m[i][j]=in[i][j];
  }
  return *this;
} 

template<int _rank>
Matrix<_rank> Matrix<_rank>::operator*(const double scal) 
{
  Matrix<_rank> out;
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) {
      out[i][j]=scal*p_m[i][j];
    }
  }
  return out;
}

template<int _rank>
Matrix<_rank> Matrix<_rank>::operator*(const Matrix<_rank>& in) 
{
  Matrix<_rank> out;
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) {
      out[i][j] = 0.;
      for(short int k=0; k<_rank; k++) out[i][j] += p_m[i][k]*in[k][j];
    }
  }
  return out;
}

template<int _rank>
void Matrix<_rank>::MatrixOut() const 
{
  double temp=0.;
  short int range=0, prcsn=0;
  short int io=msg->Out().precision(9);
  
  msg_Out()<<std::setiosflags(std::ios::fixed);
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) {
      if(temp<p_m[i][j]) temp=p_m[i][j];
    }
  }
  do { temp/=10.0; range+=1; } 
  while (temp>=1.0);
  
  msg_Out()/*<<double(range)<<msg_Out().precision()*/<<endl;
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<(range+12); j++) msg_Out()<<"-";
  }
  msg_Out()<<"-"<<endl;
  
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) {
      prcsn=-1;
      temp=p_m[i][j]-int(p_m[i][j]);
      temp=fabs(temp)/10.0;
      do {temp*=10; temp+=1.0e-14; temp=temp-int(temp); prcsn+=1;} 
      while((temp>1.0e-10) && prcsn<9);
      msg_Out()<<std::setw(range+prcsn+3)<<std::setprecision(prcsn);
      if(-1.0e-11<p_m[i][j] && p_m[i][j]<1.0e-11) msg_Out()<<double(0.0);
      else msg_Out()<<p_m[i][j];
      for(short int k=0; k<(9-prcsn); k++) msg_Out()<<" ";
      //msg_Out()<<std::setw(range+21)<<std::setprecision(18)<<p_m[i][j];
    }
    msg_Out()<<endl;
  }
  
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<(range+12); j++) msg_Out()<<"-";
  }
  msg_Out()<<"-"<<endl;
  
  msg_Out()<<endl;
  
  msg_Out()<<std::resetiosflags(std::ios::fixed); 
  msg->Out().precision(io);
}   

template<int _rank>
void Matrix<_rank>::NumRecipesNotation() 
{
  for (short int i=0;i<_rank;i++) p_m[i]--;
  p_m--;
}

template<int _rank>
void Matrix<_rank>::AmegicNotation() {
  p_m++;
  for (short int i=0;i<_rank;i++) p_m[i]++;
}

template<int _rank>
void Matrix<_rank>::Diagonalize(double* evalues,Matrix<_rank>& evectors) const
{
  double trace = 0.;
  int hit = 0;
  for (short int i=0;i<_rank;i++) trace += p_m[i][i];
  for (short int i=0;i<_rank;i++) {
    for (short int j=0;j<_rank;j++) {
      if (!IsZero(p_m[i][j]/trace)) {
	hit = 1;break;
      }
    }
  }
  if (hit==0) {
    for (short int i=0;i<_rank;i++) {
      evalues[i] = p_m[i][i];
      for (short int j=0;j<_rank;j++) evectors[i][j] = 0.;
      evectors[i][i] = 1.;
    }
    return;
  } 

  Matrix<_rank> dummy(*this);
  //minus 1  
  evectors.NumRecipesNotation();
  evalues--;
  
  int rot;
  
  dummy.NumRecipesNotation();
  dummy.Jacobi(evalues,evectors,&rot);
  dummy.AmegicNotation();
  
  evalues++;
  evectors.AmegicNotation();
}

template<int _rank>
void Matrix<_rank>::DiagonalizeSort(double* evalues,Matrix<_rank>& evectors) const
{
  int flips[_rank];
  Diagonalize(evalues, evectors);
  Matrix<_rank> flippit;
  Matrix<_rank> Mat_help;
  double help;
  int store;
  for (short int i=0;i<_rank;++i) flips[i]=i;
  for (short int i=0;i<_rank-1;++i) {
    for (short int j=i;j<_rank;++j) {
      if (dabs(evalues[i]) > dabs(evalues[j])) {
	help       = evalues[i];
	evalues[i] = evalues[j];
	evalues[j] = help;
	store      = flips[i];
	flips[i]   = flips[j];
	flips[j]   = store;
      }
    }
  }
  for (short int i=0;i<_rank;++i) flippit[flips[i]][i] = 1.;
  for (short int i=0;i<_rank;++i) {
    for (short int j=0;j<_rank;++j) {
      Mat_help[i][j] = 0;
      for(short int k=0; k<_rank; k++) Mat_help[i][j] += evectors[i][k]*flippit[k][j];
    }
  }

  evectors = Mat_help;
}

#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

template<int _rank>
void Matrix<_rank>::Jacobi(double d[], Matrix<_rank>& v, int *nrot) const
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c;
  
  double *b = new double[_rank+1];
  double *z = new double[_rank+1];
  for (ip=1;ip<=_rank;ip++) {
    for (iq=1;iq<=_rank;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=_rank;ip++) {
    b[ip]=d[ip]=p_m[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=_rank-1;ip++) {
      for (iq=ip+1;iq<=_rank;iq++)
	sm += fabs(p_m[ip][iq]);
    }
    if (sm == 0.0) {
      delete[] z;
      delete[] b;
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(_rank*_rank);
    else
      tresh=0.0;
    for (ip=1;ip<=_rank-1;ip++) {
      for (iq=ip+1;iq<=_rank;iq++) {
	g=100.0*fabs(p_m[ip][iq]);
	if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
	    && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
	  p_m[ip][iq]=0.0;
	else if (fabs(p_m[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((double)(fabs(h)+g) == (double)fabs(h))
	    t=(p_m[ip][iq])/h;
	  else {
	    theta=0.5*h/(p_m[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*p_m[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  p_m[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    ROTATE(p_m,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(p_m,ip,j,j,iq)
	      }
	  for (j=iq+1;j<=_rank;j++) {
	    ROTATE(p_m,ip,j,iq,j)
	      }
	  for (j=1;j<=_rank;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	  ++(*nrot);
	}
      }
    }
    for (ip=1;ip<=_rank;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  msg_Error()<<"Too many iterations in routine jacobi"<<endl;
}
#undef ROTATE
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */
  
template<int _rank>
Matrix<_rank> Matrix<_rank>::Dagger() 
{
  Matrix<_rank> Dag;
  for (short int i=0;i<_rank;i++) 
    for (short int j=0;j<_rank;j++) 
      Dag[i][j] = p_m[j][i];
    
  return Dag; 
}




CMatrix::CMatrix(int _rank) : m_rank(_rank) {
  p_m = new Complex*[m_rank];
  for (int i=0;i<m_rank;++i) p_m[i] = new Complex[m_rank];
  for (int i=0;i<m_rank;++i) {
    for (int j=i;j<m_rank;++j) p_m[i][j] = p_m[j][i] = Complex(0.,0.);
  }
}

CMatrix::CMatrix(Complex** m, int rank) : p_m(m), m_rank(rank) {}

CMatrix::CMatrix(const CMatrix& _in) {
  m_rank = _in.m_rank;
  p_m = new Complex*[m_rank];
  for (int i=0;i<m_rank;++i) p_m[i] = new Complex[m_rank];
  for (int i=0;i<m_rank;++i) {
    for (int j=0;j<m_rank;++j) p_m[i][j] = _in.p_m[i][j];
  }
}

CMatrix::~CMatrix() {
  for(short int i=0; i<m_rank; i++) delete[] p_m[i];
  delete[] p_m;
}


CMatrix CMatrix::Conjugate(){
  CMatrix conju = CMatrix(m_rank);
  for(int i=0; i<m_rank; i++) {
     for(int j=0; j<m_rank; j++) {
       conju[i][j] = conj(p_m[j][i]);
     }
  }
  return conju;
}

Vec4C CMatrix::operator* (const Vec4C& cvec) {
  return Vec4C(
    p_m[0][0]*cvec[0] - p_m[0][1]*cvec[1] - p_m[0][2]*cvec[2] - p_m[0][3]*cvec[3],
    p_m[1][0]*cvec[0] - p_m[1][1]*cvec[1] - p_m[1][2]*cvec[2] - p_m[1][3]*cvec[3],
    p_m[2][0]*cvec[0] - p_m[2][1]*cvec[1] - p_m[2][2]*cvec[2] - p_m[2][3]*cvec[3],
    p_m[3][0]*cvec[0] - p_m[3][1]*cvec[1] - p_m[3][2]*cvec[2] - p_m[3][3]*cvec[3]);
}



//=============================
//  Explicite instantiations.
//=============================
 
 
template class Matrix<2>; 
template class Matrix<3>; 
template class Matrix<4>; 
template class Matrix<5>;
template class Matrix<6>; 
