#ifndef Final_Selector_H
#define Final_Selector_H

#include "Primitive_Observable_Base.H"
#include <vector>

namespace ATOOLS {
  class Kt_Algorithm;

  struct Final_Selector_Data {
    bool   keep;
    int    min_n, max_n;
    double eta_min, eta_max;
    double et_min;
    double pt_min;
    double rmin;
    Final_Selector_Data(): keep(true),min_n(-1),max_n(-1),eta_min(0.),eta_max(0.),et_min(0.),pt_min(0.),rmin(0.) {}
    Final_Selector_Data(double _eta_min, double _eta_max, double _et_min, double _pt_min, double _rmin): 
      keep(true),min_n(-1),max_n(-1),eta_min(_eta_min),eta_max(_eta_max),et_min(_et_min),pt_min(_pt_min),rmin(_rmin) {}
  };

  typedef std::map<Flavour,Final_Selector_Data> Final_Data_Map;
  typedef std::map<std::pair<Flavour,Flavour>,Final_Selector_Data> Final_Correlator_Map;

  class Final_Selector : public Primitive_Observable_Base {
    std::string m_ilist, m_olist;
    Final_Data_Map m_fmap;
    Final_Correlator_Map m_cmap;

    Kt_Algorithm * p_ktalg;

    bool m_ownlist,m_extract;

    bool   PtSelect(const Vec4D & mom, double ptmin);
    bool   EtSelect(const Vec4D & mom, double etmin); 
    bool   EtaSelect(const Vec4D & mom, double etamin,double etamax);
    void   Select(Particle_List * jets,Final_Data_Map::iterator it);
    double DeltaR(const Vec4D & p1,const Vec4D & p2); 
    void   Select2(Particle_List * jets,Final_Correlator_Map::iterator it);
    void   SelectN(Particle_List * jets,Final_Data_Map::iterator it);   
    void   Extract(Particle_List * jets); 
 public:
    Final_Selector(const std::string & ilist,const std::string & olist);

    void AddSelector(const Flavour &, const Final_Selector_Data &);
    void AddSelector(const Flavour & ,const Flavour &, const Final_Selector_Data &);
    void AddSelector(const Flavour &, int, int);
    void AddKeepFlavour(const Flavour &);

    void Evaluate(const Blob_List &,double value, int ncount);
    Primitive_Observable_Base * Copy() const;    
    ~Final_Selector();
  };
}

#endif
