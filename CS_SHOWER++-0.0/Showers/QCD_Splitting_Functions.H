#ifndef QCD_Splitting_Functions_h
#define QCD_Splitting_Functions_h

#include "Splitting_Function_Base.H"
#include "PDF_Base.H"
#include "Flavour.H"

namespace CS_SHOWER {
  struct splittype {
    enum code {
      FF   = 11,
      FI   = 12,
      IF   = 21,
      II   = 22,
      none = 0
    };
  };

  //the final state emitter -- final state spectator splittings 
  class q_qg_FF : public Splitting_Function_Base {
    double m_Jmax;
    double J(const double,const double=0.,const double=0.);
  public:
    q_qg_FF(const ATOOLS::Flavour &);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };

  class q_gq_FF : public Splitting_Function_Base {
    double m_Jmax;
    double J(const double,const double=0.,const double=0.);
  public:
    q_gq_FF(const ATOOLS::Flavour &);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };

  class g_gg_FF : public Splitting_Function_Base {
    double m_Jmax;
    double J(const double,const double=0.,const double=0.);
  public:
    g_gg_FF();
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };

  class g_qq_FF : public Splitting_Function_Base {
    double m_Jmax;
    double J(const double,const double=0.,const double=0.);
  public:
    g_qq_FF(const ATOOLS::Flavour &);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };

  //the final state emitter -- initial state spectator splittings 
  class q_qg_FI : public Splitting_Function_Base {
    double m_Jmax;
    PDF::PDF_Base * p_pdf;
    double J(const double,const double,const double);
  public:
    q_qg_FI(const ATOOLS::Flavour &,PDF::PDF_Base*);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };
  
  class q_gq_FI : public Splitting_Function_Base {
    double m_Jmax;
    PDF::PDF_Base * p_pdf;
    double J(const double,const double,const double);
  public:
    q_gq_FI(const ATOOLS::Flavour &,PDF::PDF_Base*);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };
  
  class g_gg_FI : public Splitting_Function_Base {
    double m_Jmax;
    PDF::PDF_Base * p_pdf;
    double J(const double,const double,const double);
  public:
    g_gg_FI(PDF::PDF_Base*);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    
    double Z();
  };

  class g_qq_FI : public Splitting_Function_Base {
    double m_Jmax;
    PDF::PDF_Base * p_pdf;
    double J(const double,const double,const double);
  public:
    g_qq_FI(const ATOOLS::Flavour &,PDF::PDF_Base*);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };
  
  //the initial state emitter -- initial state spectator splittings 
  class q_qg_IF : public Splitting_Function_Base {
    double m_Jmax;
    PDF::PDF_Base * p_pdf;
    double J(const double,const double,const double);
  public:
    q_qg_IF(const ATOOLS::Flavour &,PDF::PDF_Base*);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
   
    double Z();
  };
  
  class q_gq_IF : public Splitting_Function_Base {
    double m_Jmax;
    PDF::PDF_Base * p_pdf;
    double J(const double,const double,const double);
  public:
    q_gq_IF(const ATOOLS::Flavour &,PDF::PDF_Base*);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };
  
  class g_gg_IF : public Splitting_Function_Base {
    double m_Jmax;
    PDF::PDF_Base * p_pdf;
    double J(const double,const double,const double);
  public:
    g_gg_IF(PDF::PDF_Base*);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };

  class g_qq_IF : public Splitting_Function_Base {
    double m_Jmax;
    PDF::PDF_Base * p_pdf;
    double J(const double,const double,const double);
  public:
    g_qq_IF(const ATOOLS::Flavour &,PDF::PDF_Base*);
    double operator() (const double,const double,
		       const double,const double);
    double OverIntegrated(const double,const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double,
			   const double,const double);
    double Z();
  };
}

#endif
