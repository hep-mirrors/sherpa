#ifndef Shower_Observables_H
#define Shower_Observables_H

#include "Primitive_Observable_Base.H"

#include "Jet_Finder.H"
#include "Kt_Algorithm.H"

namespace ATOOLS {

  class Jetrates;
  class Multiplicity;
  class PT_Distribution;

  struct Event_Obi {
    double              sum;
    Flavour             flav;
    int                 jet_ini;
    Jetrates          * jetrates;
    Multiplicity      * multi;
    PT_Distribution   * wz_pt;
    PT_Distribution   * jet_pt;
  };

  class Shower_Observables : public  Primitive_Observable_Base {
    int njet, nflavs, dohad;
    Event_Obi                all_obs;
    std::vector<Event_Obi>   jet_obs;
    std::vector<Event_Obi>   fl_obs;
    std::vector<Event_Obi>   fl_jet_obs;
    void InitObservables();
  public:
    Shower_Observables(int _type,double _xmin,double _xmax,int _nbins,
		      Selector_Base * _sel,int _njet=0,int _nflavs=0,int _dohad=0);
    void Evaluate(const Blob_List &,double,int);

    void EndEvaluation();
    void Output(const std::string & pname);

    void SetAnalysis(Primitive_Analysis*);

    Primitive_Observable_Base * Copy() const {
      std::cerr<<" ERROR in Shower_Observables::Copy() : not splittable "<<std::endl;
      return 0;
    }
  };

  struct Jetrate_Data {
    std::vector<int> jets; 
    std::vector<double> ys;    
  };

  std::ostream & operator<<(std::ostream &, const Jetrate_Data *);


  class Jetrates : public Primitive_Observable_Base {  
    Jet_Finder  * p_jfind;  //(2.,1);
    Kt_Algorithm * p_ktalg;
    std::vector<int> m_jets; 
    std::vector<double> m_ys;
    std::vector<Histogram*> m_histos;
    std::vector<Histogram*> m_rates;
    Jetrates      * p_partner;
    
  public:
    double m_ymin, m_ymax, m_sprime;

    Jetrates(int _type,double _xmin,double _xmax,int _nbins,
		      Selector_Base * _sel, const std::string & = std::string(""));
    Jetrates(Jetrates * _partner, std::string _prefix);
    Primitive_Observable_Base * Copy() const;

    void Evaluate(const Blob_List &,double, int);
    void Evaluate(const Particle_List & pl,double, int);

    void EndEvaluation();
    void EndEvaluation(double);

    void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);
    void Reset();
    ~Jetrates();
  };


  class Multiplicity : public Primitive_Observable_Base {  
    int            m_mode;
    Multiplicity * p_partner;
  public:
    Multiplicity(int _type,double _xmin,double _xmax,int _nbins, int mode, const std::string &);
    Multiplicity(Multiplicity * old);
    Multiplicity(Multiplicity * _partner, std::string _prefix);

    void Evaluate(const Particle_List & pl,double, int);
    void Evaluate(const Blob_List &,double, int);

    Primitive_Observable_Base * Copy() const;
  };


  struct ME_Data {
    int                   jets;
    double                sum;
    Flavour               flavs[3];
    std::string           name;
    Histogram * histo;
  };


  class ME_Rate : public Primitive_Observable_Base {  
    std::vector<ME_Data> m_all_rates;
    double m_sum;
  public:
    ME_Rate(int _type,double _xmin,double _xmax,int _nbins,std::string prefix);
    ME_Rate(ME_Rate * old);
    void Evaluate(const Blob_List & bl,double, int);

    void EndEvaluation();
    void EndEvaluation(double);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base * Copy() const {
      return new ME_Rate(m_type,m_xmin,m_xmax,m_nbins,"");
    }
    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);
    void Reset();

  };


  class PT_Distribution : public Primitive_Observable_Base {  
    std::vector<Histogram*>  m_histos;
    Flavour     m_checkfl;
    int         m_minn, m_maxn;
  public:
    PT_Distribution(int _type,double _xmin,double _xmax,int _nbins,
		    int _maxn, Flavour _fl,std::string =std::string("FinalState"));
    PT_Distribution(int _type,double _xmin,double _xmax,int _nbins,
		    int _minn,int _maxn, Flavour _fl,std::string =std::string("FinalState"));
    PT_Distribution(PT_Distribution * _partner, std::string _prefix);


    void Evaluate(const Particle_List &,double, int);
    void Evaluate(const Blob_List & blobs,double value, int);

    void EndEvaluation();
    void EndEvaluation(double);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base * Copy() const;
    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);
    void Reset();

  };

}

#endif
