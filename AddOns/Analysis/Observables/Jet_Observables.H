#ifndef Analysis_Observables_Jet_Observables_H
#define Analysis_Observables_Jet_Observables_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {




  /*
    mode: 1) min   m_minn jets, fill 1 (leading) to m_maxn jet 
          2) exact m_minn jets, fill 1 (leading) to m_maxn jet
  */

  class Jet_Observable_Base : public Primitive_Observable_Base {
  protected:
    unsigned int  m_mode;
    unsigned int  m_minn, m_maxn;
    std::vector<ATOOLS::Histogram*>  m_histos;
  public:

    Jet_Observable_Base(unsigned int type,double xmin,double xmax,int nbins,
			unsigned int mode, unsigned int minn, unsigned int maxn, 
			const std::string & listname=std::string("Jets"));
    virtual ~Jet_Observable_Base();

    void Evaluate(const ATOOLS::Particle_List & pl,double weight, double ncount);
    void Evaluate(const ATOOLS::Blob_List & blobs,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();

    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base & ob);
    void Reset();

    virtual double Calc(const ATOOLS::Particle * p) =0;
  };

  class Two_Jet_Observable_Base : public Primitive_Observable_Base {
  protected:
    unsigned int  m_mode;
    unsigned int  m_minn, m_maxn;
    double        * p_minpts, * p_maxpts;
    std::vector<ATOOLS::Histogram*>  m_histos;
  public:
    Two_Jet_Observable_Base(unsigned int type,double xmin,double xmax,int nbins,
			    unsigned int mode, unsigned int minn, unsigned int maxn, 
			    const std::string & =std::string("Jets"));
    ~Two_Jet_Observable_Base();

    void Evaluate(const ATOOLS::Particle_List &,double, double);
    void Evaluate(const ATOOLS::Blob_List & blobs,double value, double);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    
    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);
    void Reset();

    virtual double Calc(const ATOOLS::Particle *,const ATOOLS::Particle *,const int,const int) = 0;
    virtual void   SetPTRange(const unsigned int jetno,const double minpt,const double maxpt=1.e12);
  };

  class Three_Jet_Observable_Base : public Primitive_Observable_Base {
  protected:
    unsigned int  m_mode;
    unsigned int  m_minn, m_maxn;
    std::vector<ATOOLS::Histogram*>  m_histos;
  public:
    Three_Jet_Observable_Base(unsigned int type,double xmin,double xmax,int nbins,
			    unsigned int mode, unsigned int minn, unsigned int maxn, 
			    const std::string & =std::string("Jets"));
    ~Three_Jet_Observable_Base();

    void Evaluate(const ATOOLS::Particle_List &,double, double);
    void Evaluate(const ATOOLS::Blob_List & blobs,double value, double);
    
    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);
    void Reset();

    virtual double Calc(const ATOOLS::Particle *,const ATOOLS::Particle *,const ATOOLS::Particle *,
			const int,const int,const int) = 0;
  };

  class Jet_Eta_Distribution : public Jet_Observable_Base {

  public:
    Jet_Eta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			 unsigned int mode,unsigned int minn,unsigned int maxn, 
			 const std::string & =std::string("Jets"));

    double Calc(const ATOOLS::Particle * p);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_Phi_Distribution : public Jet_Observable_Base {
    
  public:
    Jet_Phi_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			 unsigned int mode,unsigned int minn,unsigned int maxn, 
			 const std::string & =std::string("Jets"));
    
    double Calc(const ATOOLS::Particle * p);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_Rapidity_Distribution : public Jet_Observable_Base {

  public:
    Jet_Rapidity_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			 unsigned int mode,unsigned int minn,unsigned int maxn, 
			 const std::string & =std::string("Jets"));

    double Calc(const ATOOLS::Particle * p);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_PT_Distribution : public Jet_Observable_Base {

  public:
    Jet_PT_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			unsigned int mode,unsigned int minn,unsigned int maxn, 
			const std::string & =std::string("Jets"));
    double Calc(const ATOOLS::Particle * p);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_IPT2_Distribution : public Jet_Observable_Base {

  public:
    Jet_IPT2_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			  unsigned int mode,unsigned int minn,unsigned int maxn, 
			  const std::string & =std::string("Jets"));
    double Calc(const ATOOLS::Particle * p);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_ET_Distribution : public Jet_Observable_Base {

  public:
    Jet_ET_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			unsigned int mode,unsigned int minn,unsigned int maxn, 
			const std::string & =std::string("Jets"));
    double Calc(const ATOOLS::Particle * p);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_E_Distribution : public Jet_Observable_Base {

  public:
    Jet_E_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			unsigned int mode,unsigned int minn,unsigned int maxn, 
			const std::string & =std::string("Jets"));
    double Calc(const ATOOLS::Particle * p);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_Differential_Rates : public Jet_Observable_Base {
    std::string m_reflistname;
  public:
    Jet_Differential_Rates(unsigned int type,double xmin,double xmax,int nbins,
			   unsigned int mode, unsigned int minn,unsigned int maxn, 
			   const std::string & =std::string("FinalState"),
			   const std::string & =std::string(""));
    double Calc(const ATOOLS::Particle * p);
    void Evaluate(const ATOOLS::Blob_List & blobs,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);

    Primitive_Observable_Base * Copy() const ;
  };

  //////////////////////////////////////////////////////////////////////



  class Jet_DeltaR_Distribution : public Two_Jet_Observable_Base {
  protected:
  public:
    Jet_DeltaR_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			    unsigned int mode, unsigned int minn, unsigned int maxn, 
			    const std::string & =std::string("Jets"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,
		const int jet1,const int jet2);
  };

  class Jet_DeltaEta_Distribution : public Two_Jet_Observable_Base {
  protected:
  public:
    Jet_DeltaEta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			      unsigned int mode, unsigned int minn, unsigned int maxn, 
			      const std::string & =std::string("Jets"));
    
    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,
		const int jet1,const int jet2);
  };

  class Jet_DeltaY_Distribution : public Two_Jet_Observable_Base {
  protected:
  public:
    Jet_DeltaY_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			    unsigned int mode, unsigned int minn, unsigned int maxn, 
			    const std::string & =std::string("Jets"));
    
    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,
		const int jet1,const int jet2);
  };

  class Jet_DeltaPhi_Distribution : public Two_Jet_Observable_Base {
  protected:
  public:
    Jet_DeltaPhi_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			      unsigned int mode, unsigned int minn, unsigned int maxn, 
			      const std::string & =std::string("Jets"));
    
    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,
		const int jet1,const int jet2);
  };

  class Jet_DiMass_Distribution : public Two_Jet_Observable_Base {
  protected:
  public:
    Jet_DiMass_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			    unsigned int mode, unsigned int minn, unsigned int maxn, 
			    const std::string & =std::string("Jets"));
    
    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,
		const int jet1,const int jet2);
  };

  class Eta_3_Prime : public Three_Jet_Observable_Base {
  public:
    Eta_3_Prime(unsigned int type,double xmin,double xmax,int nbins,
		unsigned int mode, unsigned int minn, unsigned int maxn, 
		const std::string & =std::string("Jets"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,const ATOOLS::Particle * p3,
		const int jet1,const int jet2,const int jet3);
  };

  class Y_3_Prime : public Three_Jet_Observable_Base {
  public:
    Y_3_Prime(unsigned int type,double xmin,double xmax,int nbins,
		unsigned int mode, unsigned int minn, unsigned int maxn, 
		const std::string & =std::string("Jets"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,const ATOOLS::Particle * p3,
		const int jet1,const int jet2,const int jet3);
  };

  class Phi_3_Prime : public Three_Jet_Observable_Base {
  public:
    Phi_3_Prime(unsigned int type,double xmin,double xmax,int nbins,
		unsigned int mode, unsigned int minn, unsigned int maxn, 
		const std::string & =std::string("Jets"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,const ATOOLS::Particle * p3,
		const int jet1,const int jet2,const int jet3);
  };

  class CosPhi_3_Prime : public Three_Jet_Observable_Base {
  public:
    CosPhi_3_Prime(unsigned int type,double xmin,double xmax,int nbins,
		unsigned int mode, unsigned int minn, unsigned int maxn, 
		const std::string & =std::string("Jets"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const ATOOLS::Particle * p1,const ATOOLS::Particle * p2,const ATOOLS::Particle * p3,
		const int jet1,const int jet2,const int jet3);
  };


}
#endif
