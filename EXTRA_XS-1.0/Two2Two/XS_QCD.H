#ifndef XS_QCD_H
#define XS_QCD_H

#include "Single_XS.H"

namespace EXTRAXS {

  struct sgq { // simple gauge qualifier
    enum code {
      none     = 0,
      pure_qcd = 1,
      pure_ew  = 2,
      mixed    = 3  };
  };

  class XS_pp_ffbar : public Single_XS { 
  public:

    XS_pp_ffbar(const size_t nin,const size_t nout,
		const ATOOLS::Flavour *fl);

    double operator()(double,double,double);
    bool SetColours(double s,double t,double u);

  };

  template <sgq::code mode>
  class XS_q1q2_q1q2 : public Single_XS {
  private:

    int a,p;
    double Ms,Mt,Mu,aS;

  public:

    XS_q1q2_q1q2(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };

  template <sgq::code mode>
  class XS_q1qbar1_q2qbar2 : public Single_XS {
  private:

    int a,p;
    double Ms,Mt,Mu,aS;

  public:

    XS_q1qbar1_q2qbar2(const size_t nin,const size_t nout, 
		       const ATOOLS::Flavour *fl);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };

  template <sgq::code mode>
  class XS_q1q1_q1q1 : public Single_XS {
  private:

    int a,p;
    double Ms,Mt,Mu,aS;

  public:

    XS_q1q1_q1q1(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl) ;

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };

  template <sgq::code mode>
  class XS_q1qbar1_q1qbar1 : public Single_XS {
  private:

    int a,p,r;
    double Ms,Mt,Mu,aS;

  public:

    XS_q1qbar1_q1qbar1(const size_t nin,const size_t nout, 
		       const ATOOLS::Flavour *fl) ;

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };

  class XS_q1qbar1_gg : public Single_XS {
  private:

    int a,p;
    double Ms,Mt,Mu,aS;

  public:

    XS_q1qbar1_gg(const size_t nin,const size_t nout, 
		  const ATOOLS::Flavour *fl);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };
  
  class XS_gg_q1qbar1 : public Single_XS {
  private:

    int r;
    double Ms,Mt,Mu,aS;

  public:

    XS_gg_q1qbar1(const size_t nin,const size_t nout, 
		  const ATOOLS::Flavour *fl);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };

  class XS_q1g_q1g : public Single_XS {
  private:
  
    int a,p;
    int fin_q,ini_q;
    int swap_ut;
    double Ms,Mt,Mu,aS;

  public:

    XS_q1g_q1g(const size_t nin,const size_t nout, 
	       const ATOOLS::Flavour *fl) ;

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };

  class XS_gg_gg : public Single_XS {
  private:
  
    double Ms,Mt,Mu,aS;

  public:

    XS_gg_gg(const size_t nin,const size_t nout, 
	     const ATOOLS::Flavour *fl);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };
 
}// end of namespace EXTRAXS

#endif

