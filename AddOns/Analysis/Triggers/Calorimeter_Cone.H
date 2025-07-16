#ifndef Analysis_Triggers_Calorimeter_Cone_H
#define Analysis_Triggers_Calorimeter_Cone_H

#include "AddOns/Analysis/Triggers/Primitive_Calorimeter.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Histogram.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  struct Jet_Data {
    int           i,j,orig;
    ATOOLS::Vec4D mom;
    double        et;
    Jet_Data(int _i, int _j, int _orig,const ATOOLS::Vec4D & _mom, double _et) :
      i(_i),j(_j),orig(_orig),mom(_mom),et(_et) {}
  };

  class Order_PT_JData {
  public:
    int operator()(const Jet_Data & a, const Jet_Data & b) {
      if (a.mom.PPerp2() > b.mom.PPerp2()) return 1;
      return 0;
    }
  };

  class Calorimeter_Cone : public Jet_Algorithm_Base {
    double m_dR, m_dR2, m_Etcut, m_Etstop, m_maxetajet, m_minetajet;
    int    m_neta, m_nphi, m_dneta, m_dnphi;
    double m_mineta, m_maxeta, m_delta_eta, m_delta_phi;
    int    m_etamode;

    int    ** p_jetno;

    Primitive_Calorimeter * p_calorimeter;
    
    std::vector<Jet_Data> m_jets;

    void CalcJets();
    void SortPT();
    void Test();
  public:
    Calorimeter_Cone(const double Etcut,const double etamin, 
		     const double etamax,double sep=1.);
    ~Calorimeter_Cone();

    bool ConstructJets(const ATOOLS::Particle_List * pl = 0,ATOOLS::Particle_List * jets = 0,
		       std::vector<double> * kt2 = 0,double rmin=-1.);
    void SetEtaRangeForJets(const double mineta,const double maxeta,const int mode);
    void FillShape(int jetno,ATOOLS::Histogram * histo,double weight=1.,double ncount=1);

    int    NumberOfJets();
    double Radius();
    double Et_cut();
    void SetAnalysis(Primitive_Analysis  * ana);
  };


  inline void Calorimeter_Cone::SetEtaRangeForJets(const double mineta,
						   const double maxeta,const int mode) 
  {
    if (mineta<maxeta) { m_minetajet = mineta; m_maxetajet = maxeta; }
                  else { m_minetajet = maxeta; m_maxetajet = mineta; }
    m_etamode = mode;
  }

  inline int    Calorimeter_Cone::NumberOfJets() { return m_jets.size(); }
  inline double Calorimeter_Cone::Radius()       { return m_dR; }
  inline double Calorimeter_Cone::Et_cut()       { return m_Etcut; }
  // a simplied version of the UA1  jet algorithm with jetradius and 
  //   minimum scalar transverse energy
  //        (RJET=1., EJCUT=5. FOR UA1)
  //  getjet.f obtained from MLM Feb 2004
  //    "Calorimeter simulation obtained from Frank Paige 23 March 1988"
  //    (RESEM=.15, RESHAD=.35 FOR URANIUM CALORIMETER)
  //    (RESEM=.15, RESHAD=.70 FOR IRON CALORIMETER)
  //    (RESEM=.11, RESHAD=.70 FOR CDF)

  //
}

#endif
