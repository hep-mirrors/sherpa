// mrst99 class file
// Comments to the C++ version to: jeppe.andersen@durham.ac.uk
// Original comments from the FORTRAN version:
//  C****************************************************************C
//  C								     C
//  C     This is a package for the new **corrected** MRST parton    C
//  C     distributions. The format is similar to the previous       C
//  C     (1998) MRST series.                                        C
//  C								     C
//  C     NOTE: 7 new sets are added here, corresponding to shifting C
//  C     the small x HERA data up and down by 2.5%, and by varying  C
//  C     the charm and strange distributions, and by forcing a      C
//  C     larger d/u ratio at large x.                               C
//  C								     C
//  C     As before, x times the parton distribution is returned,    C
//  C     q is the scale in GeV, MSbar factorization is assumed,     C
//  C     and Lambda(MSbar,nf=4) is given below for each set.        C
//  C								     C
//  C     NAMING SCHEME:                                             C
//  C						                     C
//  C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
//  C  ----  ---    -------             --------  -------   ------   C
//  C								     C
//  C  1     COR01  central gluon, a_s    300      0.1175   0.00524  C
//  C  2     COR02  higher gluon          300      0.1175   0.00497  C
//  C  3     COR03  lower gluon           300      0.1175   0.00398  C
//  C  4     COR04  lower a_s             229      0.1125   0.00585  C
//  C  5     COR05  higher a_s            383      0.1225   0.00384  C
//  C  6     COR06  quarks up             303.3    0.1178   0.00497  C
//  C  7     COR07  quarks down           290.3    0.1171   0.00593  C
//  C  8     COR08  strange up            300      0.1175   0.00524  C
//  C  9     COR09  strange down          300      0.1175   0.00524  C
//  C  10    C0R10  charm up              300      0.1175   0.00525  C
//  C  11    COR11  charm down            300      0.1175   0.00524  C
//  C  12    COR12  larger d/u            300      0.1175   0.00515  C
//  C						                     C
//  C      The corresponding grid files are called cor01.dat etc.    C
//  C							  	     C
//  C      The reference is:                                         C
//  C      A.D. Martin, R.G. Roberts, W.J. Stirling, R.S Thorne      C
//  C      Univ. Durham preprint DTP/99/64, hep-ph/9907231 (1999)    C
//  C                                                                C
//  C      Comments to : W.J.Stirling@durham.ac.uk                   C
//  C                                                                C
//  C								     C
//  C****************************************************************C

#include <iostream>
#include <fstream>
#include <string>

// old style (gcc 2.95)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ANSI-C++ style (gcc 3.2, sgi CC)

#include "mrst99.h"

using namespace PDF;
using namespace PDF::MRST99;
using namespace std;

int c_mrst99function::initialise(int mode,string path)
{
  // Initialising the x array common to all members of the class
  // Unfortunately, ANSI-C++ does not allow a initialisation of a 
  // member array. So here we go...

  xx[0]=0;
  xx[1]=1E-5;
  xx[2]=2E-5;
  xx[3]=4E-5;
  xx[4]=6E-5;
  xx[5]=8E-5;
  xx[6]=1E-4;
  xx[7]=2E-4;
  xx[8]=4E-4;
  xx[9]=6E-4;
  xx[10]=8E-4;
  xx[11]=1E-3;
  xx[12]=2E-3;
  xx[13]=4E-3;
  xx[14]=6E-3;
  xx[15]=8E-3;
  xx[16]=1E-2;
  xx[17]=1.4E-2;
  xx[18]=2E-2;
  xx[19]=3E-2;
  xx[20]=4E-2;
  xx[21]=6E-2;
  xx[22]=8E-2;
  xx[23]=.1E0;
  xx[24]=.125E0;
  xx[25]=.15E0;
  xx[26]=.175E0;
  xx[27]=.2E0;
  xx[28]=.225E0;
  xx[29]=.25E0;
  xx[30]=.275E0;
  xx[31]=.3E0;
  xx[32]=.325E0;
  xx[33]=.35E0;
  xx[34]=.375E0;
  xx[35]=.4E0;
  xx[36]=.425E0;
  xx[37]=.45E0;
  xx[38]=.475E0;
  xx[39]=.5E0;
  xx[40]=.525E0;
  xx[41]=.55E0;
  xx[42]=.575E0;
  xx[43]=.6E0;
  xx[44]=.65E0;
  xx[45]=.7E0;
  xx[46]=.75E0;
  xx[47]=.8E0;
  xx[48]=.9E0;
  xx[49]=1E0;
  
  // ditto for qq array
  
  qq[0]=0;
  qq[1]=1.25;
  qq[2]=1.5E0;
  qq[3]=2E0;
  qq[4]=2.5E0;
  qq[5]=3.2E0;
  qq[6]=4E0;
  qq[7]=5E0;
  qq[8]=6.4E0;
  qq[9]=8E0;
  qq[10]=1E1;
  qq[11]=1.2E1;
  qq[12]=1.8E1;
  qq[13]=2.6E1;
  qq[14]=4E1;
  qq[15]=6.4E1;
  qq[16]=1E2;
  qq[17]=1.6E2;
  qq[18]=2.4E2;
  qq[19]=4E2;
  qq[20]=6.4E2;
  qq[21]=1E3;
  qq[22]=1.8E3;
  qq[23]=3.2E3;
  qq[24]=5.6E3;
  qq[25]=1E4;
  qq[26]=1.8E4;
  qq[27]=3.2E4;
  qq[28]=5.6E4;
  qq[29]=1E5;
  qq[30]=1.8E5;
  qq[31]=3.2E5;
  qq[32]=5.6E5;
  qq[33]=1E6;
  qq[34]=1.8E6;
  qq[35]=3.2E6;
  qq[36]=5.6E6;
  qq[37]=1E7;
  
  // ditto for n0 array

  n0[0]=0.0;
  n0[1]=3.0;
  n0[2]=4.0;
  n0[3]=5.0;
  n0[4]=9.0;
  n0[5]=9.0;
  n0[6]=9.0;
  n0[7]=9.0;
  n0[8]=9.0;

  // The name of the file to open is stored in 'filename'
  std::cout<<"Initialise MRST99 from "<<path<<"  ("<<path.length()<<")"<<endl;

  char filename[200];
  sprintf(filename,(path+string("/cor%02d.dat")).c_str(),mode);


  std::cout<<"Initialise MRST99 from "<<path<<"  ("<<path.length()<<")  -> "
	   <<filename<<endl;

  ifstream data_file;
  data_file.open(filename);
  
  if (data_file.bad()) {
    std::cerr << "Error opening " << filename << "\n";
    exit (-1);
  }

  int i,j,n,m,k; // counters
  for (n=1;n<=nx-1;n++) 
    for (m=1;m<=nq;m++) {
      // notation: 1=uval 2=val; 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
      data_file >> f[1][n][m];
      data_file >> f[2][n][m];
      data_file >> f[3][n][m];
      data_file >> f[4][n][m];
      data_file >> f[5][n][m];
      data_file >> f[7][n][m];
      data_file >> f[6][n][m];
      data_file >> f[8][n][m];

      for (i=1;i<=np;i++) 
	f[i][n][m]=f[i][n][m]/pow((1.0-xx[n]),n0[i]);
    }      
  // close the datafile
  data_file.close();

  for (j=1;j<=ntenth-1;j++) {
    xx[j]=log10(xx[j]/xx[ntenth])+xx[ntenth];
    for (i=1;i<=8;i++) {
      if ((i==5)||(i==7))
	// skip i==5 and i==7
	i++;
      for (k=1;k<=nq;k++) 
	f[i][j][k]=log10(f[i][j][k]/f[i][ntenth][k])+f[i][ntenth][k];
    }
  }

  // zero remaining elements

  for (i=1;i<=np;i++) {
    for (m=1;m<=nq;m++) f[i][nx][m]=0.0;
  }

  std::cout<<"Initialised MRST99 from "<<path<<"  ("<<path.length()<<")  -> "<<filename<<endl;

  return 0;
}

struct s_partoncontent c_mrst99function::update(double x, double qsq)
  // Returns the parton content at x,qsq
{
  struct s_partoncontent partcontent;
  double xxx,g[np+1];
  int i,n,m;
  double a,b,fac;

  if (x<xmin) {
    std::cout << "   WARNING:  x   VALUE IS OUT OF RANGE :";
    std::cout<<x<<" not in "<<xmin<<" ... "<<xmax<<endl;
    x=xmin;
  }
  else if (x>xmax) {
    std::cout << "   WARNING:  x   VALUE IS OUT OF RANGE :";
    std::cout<<x<<" not in "<<xmin<<" ... "<<xmax<<endl;
    x=xmax;
  }
  if (qsq<qsqmin) { 
    std::cout << "   WARNING:  Q^2 VALUE IS OUT OF RANGE :";
    std::cout<<qsq<<" not in "<<qsqmin<<" ... "<<qsqmax<<endl;
    qsq=qsqmin;
  }
  else if (qsq>qsqmax) {
    std::cout << "   WARNING:  Q^2 VALUE IS OUT OF RANGE :";
    std::cout<<qsq<<" not in "<<qsqmin<<" ... "<<qsqmax<<endl;
    qsq=qsqmax;
  }

  
  xxx=x;
  
  if (x<xx[ntenth]) xxx=log10(x/xx[ntenth])+xx[ntenth];
  
  n = 1;
  while (xxx>xx[n+1]) n++;
  a = (xxx-xx[n])/(xx[n+1]-xx[n]);
  
  m=1;
  while (qsq>qq[m+1]) m++;
  b=(qsq-qq[m])/(qq[m+1]-qq[m]);

  for (i=1;i<=np;i++) {
    g[i]=(1.0-a)*(1.0-b)*f[i][n][m] + (1.0-a)*b*f[i][n][m+1] + a*(1.0-b)*f[i][n+1][m] + a*b*f[i][n+1][m+1];
    if ((n<ntenth)&&(i!=5)&&(i!=7)) {
      fac=(1.0-b)*f[i][ntenth][m]+b*f[i][ntenth][m+1];
      g[i]=fac*pow(10.0,(g[i]-fac));
    }
    g[i]=g[i]*pow((1.0-x),n0[i]);
  }
  
  partcontent.upv  = g[1];
  partcontent.dnv  = g[2];
  partcontent.usea = g[4];
  partcontent.dsea = g[8];
  partcontent.str  = g[6];
  partcontent.chm  = g[5];
  partcontent.glu  = g[3];
  partcontent.bot  = g[7];

  return partcontent;
}

c_mrst::c_mrst(std::string path)
{
  std::cout<<"Initialise MRST99 from "<<path<<"  ("<<path.length()<<")"<<endl;
  int i; // counter
  for (i=1;i<=12;i++) function[i-1].initialise(i,path);
  table[0] = &cont.upv;
  table[1] = &cont.dnv;
  table[2] = &cont.usea;
  table[3] = &cont.dsea;
  table[4] = &cont.str;
  table[5] = &cont.chm;
  table[6] = &cont.bot;
  table[7] = &cont.glu;
}

void c_mrst::mrst99(double x,double q2,int mode)
{  
  if ((mode<1)||(mode>12)) {
    std::cout << "   WARNING:  mode   VALUE IS OUT OF RANGE\n";
    mode=1;
  }

  // To save memory, we use the first array element also
  // instead of following the FORTRAN standard
  cont = function[mode-1].update(x,q2);
}
