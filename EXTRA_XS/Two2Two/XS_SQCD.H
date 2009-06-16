#ifndef EXTRA_XS_Two2Two_XS_SQCD_H
#define EXTRA_XS_Two2Two_XS_SQCD_H

#include "EXTRA_XS/Main/Single_Process.H"

namespace EXTRAXS {

  class XS_q1q2_sQ1sQ2 : public Single_Process {
  private:

    int    m_a, m_p;
    double m_mgluino2, m_g3;
    
  public:

    XS_q1q2_sQ1sQ2(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_q1q2_sQ1LsQ2R : public Single_Process {
  private:

    int    m_a, m_p;
    double m_mgluino2, m_msq32, m_msq42, m_g3;
    
  public:

    XS_q1q2_sQ1LsQ2R(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_q1qbar2_sQ1sQbar2 : public Single_Process {
  private:

    int    m_a, m_p;
    double m_mgluino2, m_msq32, m_msq42, m_g3;
    
  public:

    XS_q1qbar2_sQ1sQbar2(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_q1qbar1_sQ2sQbar2 : public Single_Process {
  private:

    int    m_a, m_p;
    double m_msquark2, m_g3;

  public:

    XS_q1qbar1_sQ2sQbar2(const size_t nin,const size_t nout, 
		       const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_q1q1_sQ1sQ1 : public Single_Process {
  private:

    int    m_a;
    double m_msquark2, m_mgluino2, m_g3;

  public:

    XS_q1q1_sQ1sQ1(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl,XS_Model_Base *const model) ;

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_q1q1_sQ1LsQ1R : public Single_Process {
  private:

    int    m_a;
    double m_msq32, m_msq42, m_mgluino2, m_g3;

  public:

    XS_q1q1_sQ1LsQ1R(const size_t nin,const size_t nout, 
		 const ATOOLS::Flavour *fl,XS_Model_Base *const model) ;

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_q1qbar1_sQ1sQbar1 : public Single_Process {
  private:

    int    m_a, m_p, m_r;
    double m_mgluino2, m_msquark2, m_g3;

  public:

    XS_q1qbar1_sQ1sQbar1(const size_t nin,const size_t nout, 
		       const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_q1qbar1_GluinoGluino : public Single_Process {
  private:

    int    m_a, m_p;
    double m_msqL2, m_msqR2, m_mgluino2, m_g3;

  public:

    XS_q1qbar1_GluinoGluino(const size_t nin,const size_t nout, 
		  const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    bool SetColours();
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };
  
  class XS_gg_sQ1sQbar1 : public Single_Process {
  private:

    int    m_r;
    double m_msquark2, m_g3;

  public:

    XS_gg_sQ1sQbar1(const size_t nin,const size_t nout, 
		  const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_q1g_sQ1Gluino : public Single_Process {
  private:
  
    int    m_a, m_p, m_finq, m_iniq, m_swaput;
    double m_msquark2, m_mgluino2, m_g3;
    
  public:

    XS_q1g_sQ1Gluino(const size_t nin,const size_t nout, 
	       const ATOOLS::Flavour *fl,XS_Model_Base *const model) ;

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };

  class XS_gg_GluinoGluino : public Single_Process {
  private:
  
    double m_g3;
    double m_mgluino2;
  public:

    XS_gg_GluinoGluino(const size_t nin,const size_t nout, 
	     const ATOOLS::Flavour *fl,XS_Model_Base *const model);

    double operator()(double s,double t,double u);
    bool SetColours(double s, double t, double u);
    friend bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm);
  };
 
  inline bool SuperPartner(const ATOOLS::Flavour & susy,const ATOOLS::Flavour & sm) {
    if (susy==ATOOLS::Flavour(kf_Gluino) && 
	sm==ATOOLS::Flavour(kf_gluon)) return true;
    if (susy.IsSquark() && sm.IsQuark() &&
	((susy.IsAnti() && sm.IsAnti()) || (!susy.IsAnti() && !sm.IsAnti()))) {
      if (abs(sm.Kfcode()) == (abs(susy.Kfcode())-1000000) ||
	  abs(sm.Kfcode()) == (abs(susy.Kfcode())-2000000)) return true;
    }
    return false;
    
  }
  
}// end of namespace EXTRAXS

#endif

