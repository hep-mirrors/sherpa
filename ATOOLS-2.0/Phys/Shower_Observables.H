#ifndef Shower_Observables_H
#define Shower_Observables_H

#include "Primitive_Observable_Base.H"

#include "Jet_Finder.H"

namespace APHYTOOLS {

  class Jetrates;
  class Multiplicity;

  struct Event_Obi {
    double              sum;
    Flavour             flav;
    int                 jet_ini;
    Jetrates          * jetrates;
    Multiplicity      * multi;
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
		      APHYTOOLS::Selector_Base * _sel,int _njet=0,int _nflavs=0,int _dohad=0);
    void Evaluate(const APHYTOOLS::Blob_List &,double);

    void EndEvaluation();
    void Output(std::string pname);
  };



  class Jetrates : public Primitive_Observable_Base {  
    APHYTOOLS::Jet_Finder  *jfind;  //(2.,1);
    std::vector<int> jets; 
    std::vector<double> ys;
    std::vector<AMATOOLS::Histogram*> histos;
    std::vector<AMATOOLS::Histogram*> rates;
    Jetrates      *partner;
    
  public:
    double ymin,ymax, sprime;

    Jetrates(int _type,double _xmin,double _xmax,int _nbins,
		      APHYTOOLS::Selector_Base * _sel);
    Jetrates::Jetrates(Jetrates * _partner, std::string _prefix);
//     Jetrates * GetCopy() {
//       return new Jetrates(this);
//     }

    void Evaluate(const APHYTOOLS::Parton_List & pl,double);

    void EndEvaluation();
    void EndEvaluation(double);

    virtual void Output(std::string pname);

    ~Jetrates() {
      delete jfind;
    }

  };


  class Multiplicity : public Primitive_Observable_Base {  
    Multiplicity * partner;
  public:
    Multiplicity(int _type,double _xmin,double _xmax,int _nbins,
			 APHYTOOLS::Selector_Base * _sel)
    {
      type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = _sel;
      name  = std::string("multi.dat");
      histo = new AMATOOLS::Histogram(type,xmin,xmax,nbins);
    };

    Multiplicity(Multiplicity * old) {
      type  = old->Type();
      sel   = old->Sel();
      xmin  = old->Xmin();
      xmax  = old->Xmax();
      nbins = old->Nbins();
      name  = old->Name();
      histo = new AMATOOLS::Histogram(old->Histo());
      histo -> Reset();
    }

    Multiplicity::Multiplicity(Multiplicity * _partner, std::string _prefix);

    void Evaluate(const APHYTOOLS::Parton_List & pl,double);

    void Evaluate(AMATOOLS::Vec4D *,APHYTOOLS::Flavour *,double);
  };


  struct ME_Data {
    int                   jets;
    double                sum;
    Flavour               flavs[3];
    std::string           name;
    AMATOOLS::Histogram * histo;
  };


  class ME_Rate : public Primitive_Observable_Base {  
    std::vector<ME_Data> all_rates;
    double sum;
  public:
    ME_Rate(int _type,double _xmin,double _xmax,int _nbins,std::string prefix)
    {
      type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = 0;
      name  = prefix+std::string("_rates.dat");
      histo = new AMATOOLS::Histogram(type,xmin,xmax,nbins);
      sum =0.;
    };

    ME_Rate(ME_Rate * old) {
      type  = old->Type();
      sel   = old->Sel();
      xmin  = old->Xmin();
      xmax  = old->Xmax();
      nbins = old->Nbins();
      name  = old->Name();
      histo = new AMATOOLS::Histogram(old->Histo());
      histo -> Reset();
    }

    void Evaluate(const APHYTOOLS::Blob_List & bl,double);

    void EndEvaluation();
    void EndEvaluation(double);
    virtual void Output(std::string pname);

  };


}

#endif
