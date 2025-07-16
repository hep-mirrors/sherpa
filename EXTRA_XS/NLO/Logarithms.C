#include "EXTRA_XS/NLO/Logarithms.H"

using namespace EXTRAXS;
using namespace ATOOLS;

Complex EXTRAXS::L0(const double &x,const double &y)
{
  double denom=1.0-x/y;
  if (dabs(denom)<1.e-7) return -1.0-denom*(0.5+denom/3.0);
  return LnRat(x,y)/denom;
}

Complex EXTRAXS::L1(const double &x,const double &y)
{
  double denom=1.0-x/y;
  if (dabs(denom)<1.e-7) return -0.5-denom/3.0*(1.0+0.75*denom);
  return (L0(x,y)+1.0)/denom;
}

Complex EXTRAXS::L2(const double &x,const double &y)
{
  double r=x/y, denom=1.0-r;
  if (dabs(denom)<1.0e-7) return (10.0+denom*(15.0+18.0*denom))/60.0;
  return (LnRat(x,y)-(0.5*(r-1.0/r)))/pow(denom,3);
}

Complex EXTRAXS::Ls0(const double &x1,const double &y1,
		     const double &x2,const double &y2)
{
  double r1=x1/y1, r2=x2/y2;
  return Lsm1(x1,y1,x2,y2)/(1.0-r1-r2);
}

Complex EXTRAXS::Ls1(const double &x1,const double &y1,
		     const double &x2,const double &y2)
{
  double r1=x1/y1, r2=x2/y2;
  return (Ls0(x1,y1,x2,y2)+L0(x1,y1)+L0(x2,y2))/(1.0-r1-r2);
}

Complex EXTRAXS::Lsm1(const double &x1,const double &y1,
		      const double &x2,const double &y2)
{
  double r1=x1/y1, r2=x2/y2;
  double omr1=1.0-r1, omr2=1.0-r2;
  Complex dilog1, dilog2;
  if (omr1>1.0) dilog1=(sqr(M_PI)/6.0-DiLog(r1))-LnRat(x1,y1)*log(omr1);
  else dilog1=DiLog(omr1);
  if (omr2>1.0) dilog2=(sqr(M_PI)/6.0-DiLog(r2))-LnRat(x2,y2)*log(omr2);
  else dilog2=DiLog(omr2);
  return dilog1+dilog2+LnRat(x1,y1)*LnRat(x2,y2)-sqr(M_PI)/6.0;
}
