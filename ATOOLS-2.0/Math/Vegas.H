#ifndef Vegas_h
#define Vegas_h

#define NDMX 50
#define MXDIM 10

namespace AMATOOLS {

  class Vegas {
    int     i,it,j,k,mds,nd,ndo,ng,npg,ia[MXDIM+1],kg[MXDIM+1];
    double  calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo;
    double  d[NDMX+1][MXDIM+1],di[NDMX+1][MXDIM+1],dt[MXDIM+1],
            dx[MXDIM+1], r[NDMX+1],x[MXDIM+1],xi[MXDIM+1][NDMX+1],xin[NDMX+1];
    double  schi,si,swgt;

    void Rebin(double rc, int nd, double * r, double * xin, double * xi);
  public:
    void Do(double * regn, int ndim, double (*fxn)(double *, double),
	    unsigned long ncall, int itmx, double *tgral, double *sd,
	    double * chi2a); 
  };
}

#endif
