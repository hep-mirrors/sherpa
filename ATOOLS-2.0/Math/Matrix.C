#include "Matrix.H"
#include "Message.H"
#include <iomanip>


using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;

template<int _rank>
Matrix<_rank>::Matrix() 
{
  m = new double*[_rank];
  for(short int i=0; i<_rank; i++) {
    m[i] = new double[_rank];
    for(short int j=0; j<_rank; j++) 
      m[i][j]=0.0;
  }
}

template<int _rank>
Matrix<_rank>::Matrix(const double ma[_rank][_rank]) 
{
  m = new double*[_rank];
  for(short int i=0; i<_rank; i++) {
    m[i] = new double[_rank];
    for(short int j=0; j<_rank; j++) 
      m[i][j]=ma[i][j];
  }
}

template<int _rank>
Matrix<_rank>::Matrix(const Matrix<_rank>& in) 
{
  m = new double*[_rank];
  for(short int i=0; i<_rank; i++) {
    m[i] = new double[_rank];
    for(short int j=0; j<_rank; j++) 
      m[i][j]=in[i][j];
  }
}

template<int _rank>
Matrix<_rank>::~Matrix() 
{
  for(short int i=0; i<_rank; i++) delete[] m[i];
  delete[] m;
}

template<int _rank>
Matrix<_rank>& Matrix<_rank>::operator=(const Matrix<_rank>& in) 
{
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) 
      m[i][j]=in[i][j];
  }
  return *this;
} 

template<int _rank>
Matrix<_rank> Matrix<_rank>::operator*(const double scal) 
{
  Matrix<_rank> out;
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) {
      out[i][j]=scal*m[i][j];
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
      for(short int k=0; k<_rank; k++) out[i][j] += m[i][k]*in[k][j];
    }
  }
  return out;
}

template<int _rank>
void Matrix<_rank>::MatrixOut() const 
{
  double temp;
  short int range=0, prcsn=0;
  short int io=msg.Out().precision(9);
  
  msg.Out()<<std::setiosflags(std::ios::fixed);
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) {
      if(temp<m[i][j]) temp=m[i][j];
    }
  }
  do { temp/=10.0; range+=1; } 
  while (temp>=1.0);
  
  msg.Out()/*<<double(range)<<msg.Out().precision()*/<<endl;
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<(range+12); j++) msg.Out()<<"-";
  }
  msg.Out()<<"-"<<endl;
  
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<_rank; j++) {
      prcsn=-1;
      temp=m[i][j]-int(m[i][j]);
      temp=fabs(temp)/10.0;
      do {temp*=10; temp+=1.0e-14; temp=temp-int(temp); prcsn+=1;} 
      while((temp>1.0e-10) && prcsn<9);
      msg.Out()<<std::setw(range+prcsn+3)<<std::setprecision(prcsn);
      if(-1.0e-11<m[i][j] && m[i][j]<1.0e-11) msg.Out()<<double(0.0);
      else msg.Out()<<m[i][j];
      for(short int k=0; k<(9-prcsn); k++) msg.Out()<<" ";
      //msg.Out()<<std::setw(range+21)<<std::setprecision(18)<<m[i][j];
    }
    msg.Out()<<endl;
  }
  
  
  for(short int i=0; i<_rank; i++) {
    for(short int j=0; j<(range+12); j++) msg.Out()<<"-";
  }
  msg.Out()<<"-"<<endl;
  
  msg.Out()<<endl;
  
  msg.Out()<<std::resetiosflags(std::ios::fixed); 
  msg.Out().precision(io);
}   

template<int _rank>
void Matrix<_rank>::NumRecipesNotation() 
{
  for (short int i=0;i<_rank;i++) m[i]--;
  m--;
}

template<int _rank>
void Matrix<_rank>::AmegicNotation() {
  m++;
  for (short int i=0;i<_rank;i++) m[i]++;
}

template<int _rank>
void Matrix<_rank>::Diagonalize(double* evalues,Matrix<_rank>& evectors)
{
  double trace = 0.;
  int hit = 0;
  for (short int i=0;i<_rank;i++) trace += m[i][i];
  for (short int i=0;i<_rank;i++) {
    for (short int j=0;j<_rank;j++) {
      if (!IsZero(m[i][j]/trace)) {
	hit = 1;break;
      }
    }
  }
  if (hit==0) {
    for (short int i=0;i<_rank;i++) {
      evalues[i] = m[i][i];
      for (short int j=0;j<_rank;j++) evectors[i][j] = 0.;
      evectors[i][i] = 1.;
    }
    return;
  } 
  Matrix<_rank> Save(*this);
  //minus 1  
  NumRecipesNotation();
  evectors.NumRecipesNotation();
  evalues--;
  
  int rot;
  
  Jacobi(evalues,evectors,&rot);
  
  evalues++;
  AmegicNotation();
  evectors.AmegicNotation();
  for (short int i=0;i<_rank;i++)
    for (short int j=0;j<_rank;j++) m[i][j] = Save[i][j];
}

template<int _rank>
void Matrix<_rank>::DiagonalizeSort(double* evalues,Matrix<_rank>& evectors)
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
  //flippit.MatrixOut();
  //msg.Out()<<"EV's"<<endl;
  //evectors.MatrixOut();
  //Mat_help = evectors*flippit;
  for (short int i=0;i<_rank;++i) {
    for (short int j=0;j<_rank;++j) {
      Mat_help[i][j] = 0;
      for(short int k=0; k<_rank; k++) Mat_help[i][j] += evectors[i][k]*flippit[k][j];
    }
  }

  evectors = Mat_help;
  //msg.Out()<<"EV's rot"<<endl;
  //evectors.MatrixOut();
}

#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

template<int _rank>
void Matrix<_rank>::Jacobi(double d[], Matrix<_rank>& v, int *nrot)
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
    b[ip]=d[ip]=m[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=_rank-1;ip++) {
      for (iq=ip+1;iq<=_rank;iq++)
	sm += fabs(m[ip][iq]);
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
	g=100.0*fabs(m[ip][iq]);
	if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
	    && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
	  m[ip][iq]=0.0;
	else if (fabs(m[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((double)(fabs(h)+g) == (double)fabs(h))
	    t=(m[ip][iq])/h;
	  else {
	    theta=0.5*h/(m[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*m[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  m[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    ROTATE(m,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(m,ip,j,j,iq)
	      }
	  for (j=iq+1;j<=_rank;j++) {
	    ROTATE(m,ip,j,iq,j)
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
  msg.Error()<<"Too many iterations in routine jacobi"<<endl;
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
      Dag[i][j] = m[j][i];
    
  return Dag; 
}

//=============================
//  Explicite instantiations.
//=============================
 
 
template class Matrix<2>; 
template class Matrix<3>; 
template class Matrix<4>; 
template class Matrix<6>; 
