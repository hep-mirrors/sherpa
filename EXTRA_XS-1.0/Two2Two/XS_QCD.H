#ifndef XS_QCD_H
#define XS_QCD_H

#include "Single_XS.H"

namespace EXTRAXS {

  class XS_pp_ffbar : public Single_XS {  // == XS_ffbar_ee but not XS_ffbar_f'fbar' !
  public:
    XS_pp_ffbar(int _nin,int _nout,APHYTOOLS::Flavour * _fl) :
      Single_XS(_nin,_nout,_fl)  
    {
      for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
      if (p_fl[2].Strong()) {
	int a  = p_fl[2].IsAnti();
	p_colours[2][a] = 500;
	p_colours[3][1-a] = 500;
      }
    }
    double operator()(double,double,double) {};
    bool SetColours(double s,double t,double u) { 
      m_scale = s;
      return 1; 
    }
  };




  class XS_q1q2_q1q2 : public Single_XS {
    int a,p;
    double Ms,Mt,Mu,aS;
  public:
    XS_q1q2_q1q2(int _nin,int _nout, 
		 APHYTOOLS::Flavour * _fl);
    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
  };

  class XS_q1qbar1_q2qbar2 : public Single_XS {
    int a,p;
    double Ms,Mt,Mu,aS;
  public:
    XS_q1qbar1_q2qbar2(int _nin,int _nout, 
		       APHYTOOLS::Flavour * _fl);
    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
  };

  class XS_q1q1_q1q1 : public Single_XS {
    int a,p;
    double Ms,Mt,Mu,aS;
  public:
    XS_q1q1_q1q1(int _nin,int _nout, 
		 APHYTOOLS::Flavour * _fl) ;
    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
  };

  class XS_q1qbar1_q1qbar1 : public Single_XS {
    int a,p;
    double Ms,Mt,Mu,aS;
  public:
    XS_q1qbar1_q1qbar1(int _nin,int _nout, 
		       APHYTOOLS::Flavour * _fl) ;
    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
  };

  class XS_q1qbar1_gg : public Single_XS {
    int a,p;
    double Ms,Mt,Mu,aS;
  public:
    XS_q1qbar1_gg(int _nin,int _nout, 
		  APHYTOOLS::Flavour * _fl);
    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
  };
  
  class XS_gg_q1qbar1 : public Single_XS {
    double Ms,Mt,Mu,aS;
  public:
    XS_gg_q1qbar1(int _nin,int _nout, 
		  APHYTOOLS::Flavour * _fl);
    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
  };

  class XS_q1g_q1g : public Single_XS {
    int a,p;
    int fin_q,ini_q;
    int swap_ut;
    double Ms,Mt,Mu,aS;
  public:
    XS_q1g_q1g(int _nin,int _nout, 
	       APHYTOOLS::Flavour * _fl) ;
    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
  };

  class XS_gg_gg : public Single_XS {
    double Ms,Mt,Mu,aS;
  public:
    XS_gg_gg(int _nin,int _nout, 
	     APHYTOOLS::Flavour * _fl);
    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
  };
 
}
#endif

