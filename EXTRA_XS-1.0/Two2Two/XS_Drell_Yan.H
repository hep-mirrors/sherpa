#ifndef XS_Drell_Yan_H
#define XS_Drell_Yan_H

#include "Single_XS.H"

namespace EXTRAXS {
  class XS_qqbar_pg : public Single_XS {
    bool barred;
  public: 
    XS_qqbar_pg(int _nin,int _nout,ATOOLS::Flavour * _fl);
    double operator()(double s,double t,double u);
    bool   SetColours(double,double,double);
  }; 

  class XS_qg_qp : public Single_XS {
    bool barred;
  public: 
    XS_qg_qp (int _nin,int _nout,ATOOLS::Flavour * _fl);
    double operator()(double s,double t,double u);
    bool   SetColours(double,double,double);
  }; 


  class XS_ee_ffbar : public Single_XS {  // == XS_ffbar_ee but not XS_ffbar_f'fbar' !
    bool   barred;
    double qe,qf,ae,af,ve,vf,mass;
    double kappa,sin2tw,MZ2,GZ2,alpha;
    double chi1,chi2,term1,term2;
    double colfac;
    int    kswitch;
    double fac,fin;
  public:
    XS_ee_ffbar(int _nin,int _nout,ATOOLS::Flavour * _fl) ;
    double operator()(double,double,double);
    bool   SetColours(double,double,double);
    /*! 
      Kfactor to first order in alphaS for this process, is included in operator()
      if kswitch = 1, 2. The choice is related to quarks in the initial or final state,
      respectively.
    */
    double KFactor(double);
  };  // end of local class definition:
 
}
#endif






