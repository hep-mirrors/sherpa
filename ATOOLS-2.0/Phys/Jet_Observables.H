#ifndef Jet_Observables_H
#define Jet_Observables_H

#include "Primitive_Observable_Base.H"

namespace ATOOLS {




  /*
    mode: 1) min   m_minn jets, fill 1 (leading) to m_maxn jet 
          2) exact m_minn jets, fill 1 (leading) to m_maxn jet
  */

  class Jet_Observable_Base : public Primitive_Observable_Base {
  protected:
    unsigned int  m_mode;
    unsigned int  m_minn, m_maxn;
    std::vector<Histogram*>  m_histos;
  public:

    Jet_Observable_Base(unsigned int mode,double xmin,double xmax,int nbins,
		    unsigned int minn, unsigned int maxn, const std::string & =std::string("FinalState"));

    void Evaluate(const Particle_List &,double, int);
    void Evaluate(const Blob_List & blobs,double value, int);

    void EndEvaluation(double);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);
    void Reset();

    virtual double Calc(const Particle *) =0;
  };


  class Jet_Eta_Distribution : public Jet_Observable_Base {

  public:
    Jet_Eta_Distribution(unsigned int mode,double xmin,double xmax,int nbins,
		    unsigned int minn,unsigned int maxn, const std::string & =std::string("FinalState"));

    double Calc(const Particle *);
    Primitive_Observable_Base * Copy() const ;
  };

  class Jet_PT_Distribution : public Jet_Observable_Base {

  public:
    Jet_PT_Distribution(unsigned int mode,double xmin,double xmax,int nbins,
		    unsigned int minn,unsigned int maxn, const std::string & =std::string("FinalState"));
    double Calc(const Particle *);
    Primitive_Observable_Base * Copy() const ;
  };

}
#endif
