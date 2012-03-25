#ifndef Analysis_Triggers_Final_Selector_H
#define Analysis_Triggers_Final_Selector_H


#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "AddOns/Analysis/Triggers/Calorimeter_Cone.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "AddOns/Analysis/Triggers/Durham_Algorithm.H"
#include "ATOOLS/Org/Smart_Pointer.H"
#include <vector>

namespace ANALYSIS {
  struct Final_Selector_Data {
    bool   keep, ko;
    int    bf, min_n, max_n;
    double eta_min, eta_max, et_min, pt_min;
    double r_min, mass_min, mass_max, f;
    Final_Selector_Data(): 
      keep(true), ko(false), bf(true), min_n(-1), max_n(-1), 
      eta_min(0.), eta_max(0.), et_min(0.), pt_min(0.),
      r_min(0.), mass_min(-1.), mass_max(-1.), f(0.5) {}
    Final_Selector_Data(double _eta_min, double _eta_max, double _et_min, double _pt_min, 
			double _r_min, double _mass_min, double _mass_max): 
      keep(true), ko(false), bf(true), min_n(-1), max_n(-1),
      eta_min(_eta_min), eta_max(_eta_max), et_min(_et_min), pt_min(_pt_min),
      r_min(_r_min), mass_min(_mass_min), mass_max(_mass_max), f(0.5) {}
  };


  typedef std::map<ATOOLS::Flavour,Final_Selector_Data>                            Final_Data_Map;
  typedef std::map<std::pair<ATOOLS::Flavour,ATOOLS::Flavour>,Final_Selector_Data> Final_Correlator_Map;



  class Final_Selector : public Analysis_Object {
  protected:
    SP(ATOOLS::Particle_Qualifier_Base) p_qualifier;
    std::string                       m_inlistname, m_outlistname;
    bool                              m_ownlist, m_extract;
    int                               m_mode;
    Jet_Algorithm_Base              * p_jetalg;
    Final_Data_Map                    m_fmap;
    Final_Correlator_Map              m_cmap;

    bool   PtSelect(const ATOOLS::Vec4D &, double);
    bool   EtSelect(const ATOOLS::Vec4D &, double); 
    bool   EtaSelect(const ATOOLS::Vec4D &, double,double);
    bool   DeltaRSelect(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double); 
    bool   MassSelect(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double,double); 
    double DeltaR(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &); 
    void   JetSelect(Particle_List * pl,const Flavour& jf);

    void   Select(ATOOLS::Particle_List *,Final_Data_Map::iterator);
    void   Select2(ATOOLS::Particle_List *,Final_Correlator_Map::iterator);
    void   SelectN(ATOOLS::Particle_List *,Final_Data_Map::iterator);   
    void   Extract(ATOOLS::Particle_List *); 
 public:
    Final_Selector(const std::string &,const std::string &);
    Final_Selector(const std::string &,const std::string &, int mode, 
		   SP(ATOOLS::Particle_Qualifier_Base) qualifier=NULL);

    void AddSelector(const ATOOLS::Flavour &, const Final_Selector_Data &);
    void AddSelector(const ATOOLS::Flavour & ,const ATOOLS::Flavour &, const Final_Selector_Data &);
    void AddSelector(const ATOOLS::Flavour &, int, int);
    void AddSelector(const ATOOLS::Flavour &, const Final_Selector_Data &,
		     Calorimeter_Cone * const);

    void AddKeepFlavour(const ATOOLS::Flavour &);

    void Evaluate(const ATOOLS::Blob_List &,double, double);
    void Output();
    Jet_Algorithm_Base * GetJetAlgorithm() { return  p_jetalg; }

    Analysis_Object * GetCopy() const;    
    ~Final_Selector();
    void SetAnalysis(Primitive_Analysis  * ana);
  };

  class Leading_Particle : public Analysis_Object {
  private:
    ATOOLS::Particle_Qualifier_Base * p_qualifier;
    std::string                       m_inlistname, m_outlistname;
    int                               m_mode;
  public:
    Leading_Particle(const std::string &,const std::string &);
    Leading_Particle(const std::string &,const std::string &, int mode, 
		     ATOOLS::Particle_Qualifier_Base * const qualifier=NULL);
    ~Leading_Particle();
    Analysis_Object * GetCopy() const;    
    void Evaluate(const ATOOLS::Blob_List &,double, double);
    void Output();
  };
}

#endif
