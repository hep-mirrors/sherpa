#ifndef Final_Selector_H
#define Final_Selector_H


#include "Primitive_Observable_Base.H"
#include "Kt_Algorithm.H"
#include "Durham_Algorithm.H"
#include <vector>

namespace ANALYSIS {
  struct Final_Selector_Data {
    bool   keep,bf;
    int    min_n, max_n;
    double eta_min, eta_max, et_min, pt_min;
    double r_min, mass_min, mass_max;
    Final_Selector_Data(): 
      keep(true), bf(true), min_n(-1), max_n(-1), 
      eta_min(0.), eta_max(0.), et_min(0.), pt_min(0.),
      r_min(0.), mass_min(-1.), mass_max(-1.) {}
    Final_Selector_Data(double _eta_min, double _eta_max, double _et_min, double _pt_min, 
			double _r_min, double _mass_min, double _mass_max): 
      keep(true), bf(true), min_n(-1), max_n(-1),
      eta_min(_eta_min), eta_max(_eta_max), et_min(_et_min), pt_min(_pt_min),
      r_min(_r_min), mass_min(_mass_min), mass_max(_mass_max) {}
  };


  typedef std::map<ATOOLS::Flavour,Final_Selector_Data>                            Final_Data_Map;
  typedef std::map<std::pair<ATOOLS::Flavour,ATOOLS::Flavour>,Final_Selector_Data> Final_Correlator_Map;



  class Final_Selector : public Primitive_Observable_Base {
    std::string            m_inlistname, m_outlistname;

    bool m_ownlist, m_extract;
    int  m_mode;

    //    ATOOLS::Kt_Algorithm * p_ktalg;
    ATOOLS::Jet_Algorithm_Base * p_jetalg;

    Final_Data_Map         m_fmap;
    Final_Correlator_Map   m_cmap;



    bool   PtSelect(const ATOOLS::Vec4D &, double);
    bool   EtSelect(const ATOOLS::Vec4D &, double); 
    bool   EtaSelect(const ATOOLS::Vec4D &, double,double);
    bool   DeltaRSelect(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double); 
    bool   MassSelect(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double,double); 
    double DeltaR(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &); 

    void   Select(ATOOLS::Particle_List *,Final_Data_Map::iterator);
    void   Select2(ATOOLS::Particle_List *,Final_Correlator_Map::iterator);
    void   SelectN(ATOOLS::Particle_List *,Final_Data_Map::iterator);   
    void   Extract(ATOOLS::Particle_List *); 
 public:
    Final_Selector(const std::string &,const std::string &, int mode=0);

    void AddSelector(const ATOOLS::Flavour &, const Final_Selector_Data &);
    void AddSelector(const ATOOLS::Flavour & ,const ATOOLS::Flavour &, const Final_Selector_Data &);
    void AddSelector(const ATOOLS::Flavour &, int, int);
    void AddKeepFlavour(const ATOOLS::Flavour &);

    void Evaluate(const ATOOLS::Blob_List &,double, int);
    void Output();
    Primitive_Observable_Base * Copy() const;    
    ~Final_Selector();
  };
}

#endif
