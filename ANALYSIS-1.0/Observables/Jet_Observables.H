#ifndef Jet_Observables_H
#define Jet_Observables_H

#include "Primitive_Observable_Base.H"

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
			const std::string & =std::string("Jets"));

    void Evaluate(const ATOOLS::Particle_List &,double, int);
    void Evaluate(const ATOOLS::Blob_List & blobs,double value, int);

    void EndEvaluation(double);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);
    void Reset();

    virtual double Calc(const ATOOLS::Particle *) =0;
  };


  class Jet_Eta_Distribution : public Jet_Observable_Base {

  public:
    Jet_Eta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			 unsigned int mode,unsigned int minn,unsigned int maxn, 
			 const std::string & =std::string("Jets"));

    double Calc(const ATOOLS::Particle *);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_PT_Distribution : public Jet_Observable_Base {

  public:
    Jet_PT_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			unsigned int mode,unsigned int minn,unsigned int maxn, 
			const std::string & =std::string("Jets"));
    double Calc(const ATOOLS::Particle *);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_E_Distribution : public Jet_Observable_Base {

  public:
    Jet_E_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			unsigned int mode,unsigned int minn,unsigned int maxn, 
			const std::string & =std::string("Jets"));
    double Calc(const ATOOLS::Particle *);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_Differential_Rates : public Jet_Observable_Base {
  public:
    Jet_Differential_Rates(unsigned int type,double xmin,double xmax,int nbins,
			   unsigned int mode, unsigned int minn,unsigned int maxn, 
			   const std::string & =std::string("FinalState"));
    double Calc(const ATOOLS::Particle *);
    void Evaluate(const ATOOLS::Blob_List & blobs,double value, int);

    Primitive_Observable_Base * Copy() const ;
  };

}
#endif
