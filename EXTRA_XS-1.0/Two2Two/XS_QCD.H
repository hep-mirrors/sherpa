#ifndef XS_QCD_H
#define XS_QCD_H

#include "Single_XS.H"

namespace EXTRAXS {

  class XS_pp_q1qbar1 : public Single_XS {
  private:

    int     m_r;
    double  m_eq, m_m2, m_g;

  public:

    XS_pp_q1qbar1(const size_t nin,const size_t nout, 
		  const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);

  };

  class XS_q1q2_q1q2 : public Single_XS {
  private:

    int    m_a, m_p;
    double m_m12, m_m22, m_g;
    
  public:

    XS_q1q2_q1q2(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);

  };

  class XS_q1qbar1_q2qbar2 : public Single_XS {
  private:

    int    m_a, m_p;
    double m_m12, m_m32, m_g;

  public:

    XS_q1qbar1_q2qbar2(const size_t nin,const size_t nout, 
		       const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);

  };

  class XS_q1q1_q1q1 : public Single_XS {
  private:

    int    m_a;
    double m_m12, m_g;

  public:

    XS_q1q1_q1q1(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl,XS_Model_Base *const model) ;

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);

  };

  class XS_q1qbar1_q1qbar1 : public Single_XS {
  private:

    int    m_a, m_p, m_r;
    double m_m12, m_g;

  public:

    XS_q1qbar1_q1qbar1(const size_t nin,const size_t nout, 
		       const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);

  };

  class XS_q1qbar1_gg : public Single_XS {
  private:

    int    m_a, m_p;
    double m_m12, m_g;

  public:

    XS_q1qbar1_gg(const size_t nin,const size_t nout, 
		  const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();

  };
  
  class XS_gg_q1qbar1 : public Single_XS {
  private:

    int    m_r;
    double m_m32, m_g;

  public:

    XS_gg_q1qbar1(const size_t nin,const size_t nout, 
		  const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);

  };

  class XS_q1g_q1g : public Single_XS {
  private:
  
    int    m_a, m_p, m_finq, m_iniq, m_swaput;
    double m_mq2, m_g;

  public:

    XS_q1g_q1g(const size_t nin,const size_t nout, 
	       const ATOOLS::Flavour *fl,XS_Model_Base *const model) ;

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);

  };

  class XS_gg_gg : public Single_XS {
  private:
  
    double m_g;

  public:

    XS_gg_gg(const size_t nin,const size_t nout, 
	     const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);

  };
 
}// end of namespace EXTRAXS

#endif

